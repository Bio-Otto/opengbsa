"""Regression test: an unusable requested platform no longer falls back silently.

Up to 0.0.11, a requested platform that could not be used (e.g. a CUDA runtime
newer than the driver supports) fell back to the next one with only a log line,
while the startup banner still named the requested platform; a benchmark could
report a CPU timing labelled as GPU. Now the fallback is recorded
(platform_requested / platform_used / platform_fallback, also written to
results_summary.yaml), and platform_settings.strict_platform: true turns it into
an error.
"""
from types import SimpleNamespace

import openmm
import pytest

from mmgbsa.mmgbsa_core import GBSACalculator

select = GBSACalculator.setup_optimized_platform
AVAILABLE = {openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())}


@pytest.fixture(autouse=True)
def no_env_platform(monkeypatch):
    monkeypatch.delenv("OPENMM_DEFAULT_PLATFORM", raising=False)


def calc(**kw):
    base = dict(platform_preference=None, strict_platform=False)
    base.update(kw)
    return SimpleNamespace(**base)


def test_usable_request_is_used_and_recorded():
    c = calc(platform_preference="CPU")
    platform, _ = select(c)
    assert platform.getName() == "CPU"
    assert (c.platform_requested, c.platform_used, c.platform_fallback) == ("CPU", "CPU", None)


def test_unusable_request_falls_back_and_says_so(capsys):
    c = calc(platform_preference="NoSuchPlatform")
    platform, _ = select(c)
    assert platform.getName() in AVAILABLE
    assert c.platform_requested == "NoSuchPlatform"
    assert c.platform_used == platform.getName()
    assert c.platform_fallback and c.platform_fallback.startswith("NoSuchPlatform:")
    assert f"Running on OpenMM platform: {platform.getName()}" in capsys.readouterr().out


def test_strict_platform_refuses_to_fall_back():
    c = calc(platform_preference="NoSuchPlatform", strict_platform=True)
    with pytest.raises(RuntimeError, match="NoSuchPlatform"):
        select(c)


def test_auto_selection_is_recorded_as_auto():
    c = calc(platform_preference="auto")
    platform, _ = select(c)
    assert c.platform_requested == "auto" and c.platform_fallback is None
    assert c.platform_used == platform.getName()


def test_explicit_platform_name_overrides_preference():
    c = calc(platform_preference="NoSuchPlatform")
    platform, _ = select(c, platform_name="CPU")
    assert platform.getName() == "CPU" and c.platform_fallback is None
