import gerbls
import numpy as np
import numpy.typing as npt
from copy import deepcopy
from typing import Optional

# Optional dependencies
try:
    import batman
    _HAS_BATMAN = True
except ImportError:
    _HAS_BATMAN = False
try:
    from scipy.optimize import curve_fit
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False


class TransitModel:
    
    def __init__(self):

        # Parameters related to fitting
        self.chi2 = 0.
        self.chi2_const = 0.
        self.chi2r = 0.
        self.dchi2 = 0.
    
    @property
    def fitted(self):
        return (self.chi2 != 0)
    
    @staticmethod
    def get_chi2_const(phot: gerbls.pyDataContainer):
        mag0_fit = np.average(phot.mag, weights=phot.err**-2)
        return np.sum(((phot.mag - mag0_fit) / phot.err)**2)


class LDModel(TransitModel):

    def __init__(self,
                 b: float = 0.,
                 mag0: float = 1.,
                 P: float = 0.,
                 r: float = 0.,
                 t0: float = 0.,
                 target: Optional[gerbls.pyTarget] = None):
        
        super().__init__()

        self.b = b
        self.mag0 = mag0
        self.P = P
        self.r = r
        self.t0 = t0

        # Make a copy of target (if given), otherwise use a default one (with Solar values)
        self.target = (gerbls.pyTarget() if target is None else target.copy())

        # Uncertainties for the fitted parameters
        self.b_err = 0.
        self.mag0_err = 0.
        self.P_err = 0.
        self.r_err = 0.
        self.t0_err = 0.
        self.target_u1_err = 0.
        self.target_u2_err = 0.
    
    def __str__(self):
        return (
            f"{self.__class__.__name__} with the following orbital parameters:\n" +
            ("" if self.fitted else "[WARNING - parameters have not been fitted]\n") +
            f"    b = {self.b:.6f} +/- {self.b_err:.6f}\n"
            f" mag0 = {self.mag0:.6f} +/- {self.mag0_err:.6f}\n"
            f"    P = {self.P:.6f} +/- {self.P_err:.6f}\n"
            f"    r = {self.r:.6f} +/- {self.r_err:.6f}\n"
            f"   t0 = {self.t0:.6f} +/- {self.t0_err:.6f}\n"
            "and the following stellar parameters:\n"
            f"   u1 = {self.target.u1:.6f} +/- {self.target_u1_err:.6f}\n"
            f"   u2 = {self.target.u2:.6f} +/- {self.target_u2_err:.6f}\n"
            f"    M =   {self.target.M:.4f}\n"
            f"    R =   {self.target.R:.4f}\n"
            "and the following fit parameters:\n"
            f" chi2 =   {self.chi2:.4f} (reduced: {self.chi2r:.4f})\n"
            f"dchi2 =   {self.dchi2:.4f}"
        )
    
    def _get_batman_TransitParams(self):
        params = batman.TransitParams()
        params.t0 = self.t0
        params.per = self.P
        params.rp = self.r
        params.a = self.target.get_aR_ratio(self.P)
        params.inc = self.target.get_inc(self.P, self.b)
        params.ecc = 0.  # eccentricity
        params.w = 90.  # longitude of periastron (in degrees)
        params.u = self.target.u  # limb-darkening coefficients [u1, u2]
        params.limb_dark = "quadratic"  # limb darkening model
        return params

    @property
    def dur(self):
        return self.target.get_transit_duration(self.P, self.b)

    # Return a limb-darkened light curve evaluated at a given array of times
    def eval(self, time: npt.ArrayLike):

        if not _HAS_BATMAN:
            gerbls.raise_import_error("LDModel.eval", "batman")

        params = self._get_batman_TransitParams()
        return batman.TransitModel(params, np.asarray(time)).light_curve(params)

    def fit(self, phot: gerbls.pyDataContainer, u_fixed: bool = False):

        if not _HAS_BATMAN:
            gerbls.raise_import_error("LDModel.fit", "batman")
        if not _HAS_SCIPY:
            gerbls.raise_import_error("LDModel.fit", "scipy.optimize.curve_fit")
        if u_fixed and self.target.u1 == 0 and self.target.u2 == 0:
            raise ValueError(
                "Quadratic limb-darkening coefficients cannot be zero if u_fixed is True.")

        params = self._get_batman_TransitParams()
        model = batman.TransitModel(params, phot.rjd)

        # Get [u1, u2] from the [q1, q2] LD reparametrization by Kipping (2013)
        def u(q1, q2):
            return [2 * q1**0.5 * q2, q1**0.5 * (1 - 2 * q2)]
        def u_err(q1, q2, q1_err, q2_err):
            u_ = u(q1, q2)
            return [((0.5 * q1_err / q1)**2 + (q2_err / q2)**2)**0.5 * abs(u_[0]),
                    ((0.5 * u_[1] * q1_err / q1)**2 + (u_[0] * q2_err / q2)**2)**0.5]

        # Function to optimize
        # Args: [t0_phase, P, r, b, mag0] (if u_fixed)
        # Args: [t0_phase, P, r, b, mag0, u1, u2] (if not u_fixed)
        def f(x, *args):
            params.t0 = args[0] * args[1]
            params.per = args[1]
            params.rp = args[2]
            params.a = self.target.get_aR_ratio(args[1])
            params.inc = self.target.get_inc(args[1], args[3])
            if not u_fixed:
                params.u = u(*args[5:7])
            return model.light_curve(params) * args[4]

        # Initial guess, lower and upper bounds for fitted parameters
        p0 = [self.t0 / self.P, self.P, self.r, self.b, self.mag0]
        lb = [0., self.P * 0.99, 0., 0., 0.]
        ub = [1., self.P * 1.01, 1., 1., np.inf]
        if not u_fixed:
            if self.target.u1 + self.target.u2 == 0:
                p0 += [0.3, 0.3]
            else:
                u_sum = self.target.u1 + self.target.u2
                p0 += [u_sum**2, 0.5 * self.target.u1 / u_sum]
            lb += [0., 0.]
            ub += [1., 1.]

        popt, pcov = curve_fit(f, phot.rjd, phot.mag, p0=p0, sigma=phot.err, bounds=(lb, ub))
        perr = np.sqrt(np.diag(pcov))

        # Save results
        self.t0 = popt[0] * popt[1]
        self.t0_err = ((perr[0] * popt[1])**2 + (popt[0] * perr[1])**2)**0.5
        self.P = popt[1]
        self.P_err = perr[1]
        self.r = popt[2]
        self.r_err = perr[2]
        self.b = popt[3]
        self.b_err = perr[3]
        self.mag0 = popt[4]
        self.mag0_err = perr[4]
        if not u_fixed:
            self.target.u1, self.target.u2 = u(*popt[5:7])
            self.target_u1_err, self.target_u2_err = u_err(*popt[5:7], *perr[5:7])

        self.chi2 = np.sum(((f(phot.rjd, *popt) - phot.mag) / phot.err)**2)
        self.chi2_const = self.get_chi2_const(phot)
        self.chi2r = self.chi2 / (phot.size - 1)
        self.dchi2 = self.chi2 - self.chi2_const
    
    @classmethod
    def from_BLS(cls, bls: gerbls.pyBLSResult, target: Optional[gerbls.pyTarget] = None):
        b = 0. if target is None else target.estimate_b(bls.P, bls.dur)
        r = (bls.dmag / bls.mag0)**0.5
        return cls(b=b, mag0=bls.mag0, P=bls.P, r=r, t0=bls.t0, target=target)