import json

class Config:
    properties_path: str | None = None

    @staticmethod
    def selectPropertyFile(path: str) -> None:
        """Set the path to the property file."""
        Config.properties_path = path

    @staticmethod
    def get(path: str):
        """Provide the path within the properties file to the given variable."""
        with open(Config.properties_path, 'r') as file:
            config = json.load(file)
            vars = path.split('/')
            for var in vars:
                config = config[var]
            return config
