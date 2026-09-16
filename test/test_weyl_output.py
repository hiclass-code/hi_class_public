"""CLI Weyl spectra: build class and install the matching package, then run:
python -m pytest --import-mode=importlib --confcutdir=test test/test_weyl_output.py
"""
from pathlib import Path
import os
import subprocess

import numpy as np
import pytest
from hiclassy import HiClass

CLASS = Path(os.environ.get('CLASS_EXECUTABLE', Path(__file__).resolve().parents[1] / 'class'))
BASE = {'h': 0.67, 'A_s': 2.1e-9, 'n_s': 0.965, 'P_k_max_1/Mpc': 0.15}


def run_class(tmp_path, **params):
    ini = tmp_path / 'test.ini'
    config = {**BASE, 'output': 'wPk', 'root': str(tmp_path / 'result'),
              'overwrite_root': 'yes', **params}
    ini.write_text('\n'.join(f'{k} = {v}' for k, v in config.items()) + '\n')
    result = subprocess.run([str(CLASS), str(ini)], capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'Error' not in result.stdout, result.stdout
    return config


def cosmology(config):
    cosmo = HiClass()
    cosmo.set({k: v for k, v in config.items()
               if k not in {'root', 'overwrite_root', 'headers', 'format'}})
    cosmo.compute()
    return cosmo


@pytest.mark.parametrize('redshifts', ['0', '1', '0, 0.5, 1'])
def test_weyl_only_matches_python(tmp_path, redshifts):
    config = run_class(tmp_path, z_pk=redshifts)
    zs = [float(z) for z in redshifts.split(',')]
    cosmo = cosmology(config)
    try:
        for i, z in enumerate(zs):
            suffix = '' if len(zs) == 1 else f'z{i+1}_'
            path = tmp_path / f'result_{suffix}pk_weyl.dat'
            assert '2:Q_W (h/Mpc)' in path.read_text()
            data = np.loadtxt(path)
            k = data[:, 0] * cosmo.h()
            # File precision may round the last k infinitesimally above k_max.
            expected = [cosmo.pk_weyl_lin(q, z) for q in k[1:-1:5]]
            np.testing.assert_allclose(data[1:-1:5, 1] * cosmo.h(), expected, rtol=2e-10)
        assert not list(tmp_path.glob('*_pk.dat'))
    finally:
        cosmo.struct_cleanup()


def test_matter_and_weyl_against_transfer_functions(tmp_path):
    config = run_class(tmp_path, output='mPk,wPk,mTk')
    assert (tmp_path / 'result_pk.dat').exists()
    data = np.loadtxt(tmp_path / 'result_pk_weyl.dat')
    cosmo = cosmology(config)
    try:
        transfer = cosmo.get_transfer(z=0)
        k = transfer['k (h/Mpc)'] * cosmo.h()
        primordial = BASE['A_s'] * (k / 0.05) ** (BASE['n_s'] - 1)
        # Independent normalization: dimensional primordial spectrum times
        # the squared Weyl transfer, multiplied by k**4.
        expected = 2 * np.pi**2 * k * primordial * ((transfer['phi'] + transfer['psi']) / 2)**2
        np.testing.assert_allclose(data[:, 0] * cosmo.h(), k, rtol=1e-11)
        np.testing.assert_allclose(data[:, 1] * cosmo.h(), expected, rtol=2e-10)
    finally:
        cosmo.struct_cleanup()


@pytest.mark.parametrize('redshifts', ['0', '0, 1'])
def test_numbered_runs_do_not_overwrite(tmp_path, redshifts):
    for _ in range(2):
        run_class(tmp_path, z_pk=redshifts, overwrite_root='no')
    suffix = '' if redshifts == '0' else 'z1_'
    assert (tmp_path / f'result00_{suffix}pk_weyl.dat').exists()
    assert (tmp_path / f'result01_{suffix}pk_weyl.dat').exists()


def test_headers_can_be_disabled(tmp_path):
    run_class(tmp_path, headers='no', format='camb')
    path = tmp_path / 'result_pk_weyl.dat'
    assert '#' not in path.read_text()
    assert np.loadtxt(path).shape[1] == 2


def test_correlated_initial_conditions(tmp_path):
    run_class(tmp_path, ic='ad,cdi', f_cdi=0.3, c_ad_cdi=-0.4)
    load = lambda suffix: np.loadtxt(tmp_path / f'result_pk_weyl{suffix}.dat')[:, 1]
    total, ad, cdi, cross = (load(s) for s in ['', '_ad', '_cdi', '_ad_cdi'])
    np.testing.assert_allclose(total, ad + cdi + 2*cross, rtol=1e-10)
    assert np.any(cross != 0)


def test_nonlinear_matter_keeps_weyl_linear(tmp_path):
    config = run_class(tmp_path, output='mPk,wPk', **{'non linear': 'halofit', 'P_k_max_1/Mpc': 1.0})
    assert (tmp_path / 'result_pk_nl.dat').exists()
    assert not (tmp_path / 'result_pk_weyl_nl.dat').exists()
    data = np.loadtxt(tmp_path / 'result_pk_weyl.dat')
    cosmo = cosmology(config)
    try:
        indices = np.arange(1, len(data)-1, 10)
        expected = [cosmo.pk_weyl_lin(data[i, 0]*cosmo.h(), 0) for i in indices]
        np.testing.assert_allclose(data[indices, 1]*cosmo.h(), expected, rtol=2e-10)
    finally:
        cosmo.struct_cleanup()


def test_modified_gravity_weyl_output(tmp_path):
    config = run_class(tmp_path, Omega_Lambda=0, Omega_fld=0, Omega_smg=-1,
                       gravity_model='propto_omega', parameters_smg='1., 0.1, 0., 0., 1.',
                       expansion_model='lcdm', expansion_smg=0.5)
    data = np.loadtxt(tmp_path / 'result_pk_weyl.dat')
    cosmo = cosmology(config)
    try:
        expected = [cosmo.pk_weyl_lin(k*cosmo.h(), 0) for k in data[1:-1:5, 0]]
        np.testing.assert_allclose(data[1:-1:5, 1]*cosmo.h(), expected, rtol=2e-10)
    finally:
        cosmo.struct_cleanup()
