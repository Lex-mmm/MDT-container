import docker

client = docker.from_env()
container_name = 'mdt-pt-container'
image_name = 'mdt-pt'
tag = 'latest'

# Check and remove any existing container with the same name
try:
    container = client.containers.get(container_name)
    print(f"Container {container_name} already exists. Removing it...")
    container.remove(force=True)  # Remove container even if running
except docker.errors.NotFound:
    # Container doesn't exist; nothing to remove
    pass

print(f"🚀 Starting container {container_name} using image {image_name}:{tag}...")

# Now run the container safely using the proper image name
client.containers.run(
    f'{image_name}:{tag}',
    name=container_name,
    ports={'6379/tcp': 6379},
    detach=True,
    tty=True  # Optional: allocate a pseudo-TTY
)