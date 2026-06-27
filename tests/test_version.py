import packaging.version

import nrf1_proteomics


def test_version_is_valid() -> None:
    _ = packaging.version.parse(nrf1_proteomics.__version__)
