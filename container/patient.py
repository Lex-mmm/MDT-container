from randomgenerator import RandomGenerator
from config import Config
from datetime import datetime

class Patient:
    def __init__(self) -> None:
        # Generate patient specific properties
        self.id: str = str(RandomGenerator.generateRandomIdentifier(7))
        self.bsn: int = RandomGenerator.generateRandomIdentifier(9)
        self.dob: datetime = datetime.strptime(Config.get('patient/dob'), "%Y-%m-%d")
        # Additional patient-specific attributes as needed
