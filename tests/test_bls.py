import gerbls
import numpy as np

def test_bls_basic(phot_test):

    correct_P = 1.37

    results = gerbls.run_bls(phot_test.rjd, phot_test.mag, phot_test.err, 
                             min_period=0.4, max_period=10., t_samp=10/60/24)
    
    assert 'P' in results
    assert abs(results['P'][np.argmax(results['dchi2'])] - correct_P) < 0.001

def test_bls_fast(phot_test):

    correct_P = 1.37

    bls = gerbls.pyFastBLS()
    bls.create(phot_test, 0.4, 10., t_samp=10/60/24)
    bls.run()

    blsa = gerbls.pyBLSAnalyzer(bls)
    best_model = blsa.generate_models(1)[0]

    assert len(blsa.P) > 0
    assert abs(best_model.P - correct_P) < 0.001