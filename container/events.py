from datetime import datetime
import random
import logging

logging.basicConfig(level=logging.INFO)

class EventManager:
    def __init__(self, refresh_rate=1000):
        self.running = True
        self.event_pauze = False
        self.treatment_pauze = False
        self.refresh_rate = refresh_rate
        self.events = ["pneumonia", "sepsis", "ARDS"]
        self.event_info = {
            "pneumonia": {"chance": 20, "organ": "pulmonary"},
            "sepsis": {"chance": 15, "organ": "systemic"},
            "ARDS": {"chance": 10, "organ": "pulmonary"}
        }
        self.vitals = {
            "latest_epoch": {"spo2": 90, "bp": 120},
            "goal": {"spo2": {"min": 92, "max": 100}, "bp": {"min": 110, "max": 130}}
        }
        self.intervention_status = {"spo2": 0, "bp": 0}
        self.event_history = []

    def generate_events(self):
        """Generate and trigger events based on random chance."""
        if self.running and not self.event_pauze:
            for event in self.events:
                chance = self.event_info[event]['chance']
                if self.generate_single_in_range(0, 100, 10) < chance:
                    self.event_history.append({"event": event, "time": datetime.now(), "status": "untreated"})
                    logging.info(f"Event triggered: {event} at {datetime.now()}")
                    self.execute_event(event)
                    self.event_cooldown()

    def execute_event(self, event):
        """Execute event-specific actions."""
        logging.info(f"Executing event: {event} with action based on {self.event_info[event]['organ']}")

    def event_cooldown(self):
        """Cooldown for event generation to avoid rapid triggering of events."""
        self.event_pauze = True
        # Simulate cooldown (in an async or threaded app, you would reset after a delay)
        # For testing, we can simply print the cooldown status.
        logging.info("Event cooldown triggered. No events will generate until cooldown is over.")
        # Reset pause flag to allow events again (e.g., after a cooldown period in actual usage)
        self.event_pauze = False

    def generate_single_in_range(self, min_val, max_val, precision=0):
        """Generates a random integer within a given range and precision."""
        scale = 10 ** precision
        return random.randint(int(min_val * scale), int(max_val * scale)) / scale
