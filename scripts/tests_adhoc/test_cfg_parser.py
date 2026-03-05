import yaml
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mmgbsa.config import ConfigManager
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


try:
    with open("test_cfg.yaml") as f:
        y = yaml.safe_load(f)
        print("YAML RAW:", y["analysis_settings"].get("frame_start", "MISSING"))
except Exception as e:
    print("YAML LOAD ERROR:", e)

try:
    man = ConfigManager("test_cfg.yaml")
    cfg = man.get_config()
    print("CONFIG MNGR:", cfg["analysis_settings"].get("frame_start", "MISSING"))
except Exception as e:
    print("CONFIG MNGR ERROR:", e)
