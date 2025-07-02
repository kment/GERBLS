import gerbls
import numpy as np

def test_bls_basic(phot_test):
    """
    Test for the BLS basic usage function.
    """

    correct_P = 1.37

    results = gerbls.run_bls(phot_test.rjd, phot_test.mag, phot_test.err, 
                             min_period=0.4, max_period=10., t_samp=10/60/24)
    
    # Check the validity of BLS results
    assert 'P' in results
    assert len(results['P']) > 0
    assert not np.any(np.isnan(results['P']))
    assert not np.any(np.isnan(results['dchi2']))
    assert not np.any(np.isnan(results['t0']))
    assert not np.any(np.isnan(results['dur']))
    assert not np.any(np.isnan(results['mag0']))
    assert not np.any(np.isnan(results['dmag']))
    
    # Check whether the correct period was recovered
    assert abs(results['P'][np.argmax(results['dchi2'])] - correct_P) < 0.001

def test_bls_bf(phot_test):
    """
    Test for the brute-force BLS.
    """

    correct_P = 1.37

    bls = gerbls.pyBruteForceBLS()
    bls.setup(phot_test, 0.4, 10., t_bins=10/60/24)
    bls.run()

    blsa = gerbls.pyBLSAnalyzer(bls)
    assert len(blsa.P) > 0

    best_model = blsa.generate_models(1)[0]
    assert abs(best_model.P - correct_P) < 0.001

def test_bls_fast(phot_test):
    """
    Test for the fast-folded BLS.
    """

    correct_P = 1.37

    bls = gerbls.pyFastBLS()
    bls.setup(phot_test, 0.4, 10., t_samp=10/60/24)
    bls.run()

    blsa = gerbls.pyBLSAnalyzer(bls)
    assert len(blsa.P) > 0

    best_model = blsa.generate_models(1)[0]
    assert abs(best_model.P - correct_P) < 0.001