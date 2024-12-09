import random
import string
from datetime import datetime, timedelta

class RandomGenerator:
    def __init__(self, seed: str) -> None:
        """Initializes the RandomGenerator with a seed for reproducibility."""
        random.seed(seed)

    @staticmethod
    def generateRandomIdentifier(length: int) -> str:
        """Generates a random alphanumeric identifier of a given length."""
        return ''.join(random.choices(string.ascii_uppercase + string.digits, k=length))

    @staticmethod
    def generateRandomDate(start: datetime, end: datetime) -> datetime:
        """Generates a random datetime between two dates."""
        return start + timedelta(
            seconds=random.randint(0, int((end - start).total_seconds()))
        )
