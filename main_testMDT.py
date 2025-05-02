#!/usr/bin/env python3
import threading
import time
from datetime import datetime

from digital_twin_model import DigitalTwinModel

# Define a dummy redis client with a no-op add_vital_sign method
class DummyRedisClient:
    def add_vital_sign(self, patient_id, data):
        # No operation; ignore calls during test
        pass

def monitor_loop(model, run_time=500):
    """
    Let the model run for run_time seconds, printing HR and SaO2 each second.
    """
    start = time.time()
    while time.time() - start < run_time:
        time.sleep(1)
        ts = datetime.now().strftime("%H:%M:%S")

        hr  = getattr(model, "current_heart_rate", None)
        sa  = getattr(model, "current_SaO2",       None)
        ma  = getattr(model, "recent_MAP",        None)
        rr  = getattr(model, "current_RR",         None)

        # Guard against None
        hr_str = f"{hr:.1f}" if hr is not None else "n/a"
        sa_str = f"{sa:.1f}" if sa is not None else "n/a"
        ma_str = f"{ma:.1f}" if ma is not None else "n/a"
        rr_str = f"{rr:.1f}" if rr is not None else "n/a"

        print(f"[{ts}] HR={hr_str}   MAP={ma_str} RR={rr_str} SaO2={sa_str} ")

    # Stop the simulation
    model.stop_simulation()


def main():
    # 1) Instantiate your twin
    twin = DigitalTwinModel(
        patient_id="demo_patient",
        param_file="healthyFlat.json",
        sleep=True,
        time_step=0.02
    )
    # Inject the dummy redis client so that redis calls are ignored during testing
    twin.redisClient = DummyRedisClient()

    # 2) Run it in its own thread
    sim_thread = threading.Thread(target=twin.start_simulation, daemon=True)
    sim_thread.start()

    # 3) Monitor for 10 seconds then shut off
    monitor_loop(twin, run_time=500)

    sim_thread.join()
    print("Done.")


if __name__ == "__main__":
    main()