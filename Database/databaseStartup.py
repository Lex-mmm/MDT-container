## Redis in-Memory database
import docker
import time, json
import redis


### INITIALIZE TO START DOCKER CONTAINER WITH REDIS JSON CONTAINER

if __name__ == "__main__":
        
    # Initialize Docker client
    try:
        client = docker.from_env()
        print("Docker client initialized")
    except docker.errors.DockerException as e:
        print(f"Error initializing Docker client: {e}")
        exit(1)

    container_name = 'redis-stack'
    try:
        container = client.containers.get(container_name)
        print(f"RedisJSON Docker container found, status: {container.status}")
        if container.status != 'running':
            print("Starting the RedisJSON Docker container...")
            container.start()
            time.sleep(2)
            print("RedisJSON Docker container started")
        else:
            container.kill()
            container.remove()
            #container.reload()
            time.sleep(1)
            container.start()
            print("RedisJSON Docker container was already running, restarting it...")
    except docker.errors.NotFound:
        print("RedisJSON Docker container not found. Creating and starting a new one...")
        container = client.containers.run('redis/redis-stack:latest', name=container_name,
                                    ports={'6379/tcp': 6379})
        time.sleep(2)
        print("RedisJSON Docker container started")