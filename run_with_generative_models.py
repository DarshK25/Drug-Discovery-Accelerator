import os
import subprocess
import sys
import time
from threading import Thread

def install_dependencies():
    """Install the required dependencies for generative models."""
    print("Installing dependencies for generative models...")
    
    # Check if we're in a virtual environment
    in_venv = hasattr(sys, 'real_prefix') or (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix)
    
    if not in_venv:
        print("Warning: It's recommended to run this in a virtual environment.")
        response = input("Continue without a virtual environment? (y/n): ")
        if response.lower() != 'y':
            print("Exiting. Please create and activate a virtual environment first.")
            sys.exit(1)
    
    # Install backend dependencies
    os.chdir('backend')
    subprocess.run([sys.executable, '-m', 'pip', 'install', '-r', 'requirements.txt'])
    os.chdir('..')
    
    # Install frontend dependencies
    os.chdir('frontend')
    subprocess.run(['npm', 'install'])
    os.chdir('..')
    
    print("Dependencies installed successfully!")

def run_backend():
    """Run the FastAPI backend server."""
    os.chdir('backend')
    if sys.platform == 'win32':
        subprocess.run([sys.executable, '-m', 'uvicorn', 'main:app', '--reload', '--host', '0.0.0.0', '--port', '8000'])
    else:
        subprocess.run([sys.executable, '-m', 'uvicorn', 'main:app', '--reload', '--host', '0.0.0.0', '--port', '8000'])

def run_frontend():
    """Run the React frontend development server."""
    os.chdir('frontend')
    if sys.platform == 'win32':
        subprocess.run(['npm', 'start'])
    else:
        subprocess.run(['npm', 'start'])

def main():
    """Run both the frontend and backend servers with generative models."""
    print("Starting Drug Discovery Assistant with Generative Models...")
    
    # Ask if dependencies should be installed
    response = input("Do you want to install/update dependencies? (y/n): ")
    if response.lower() == 'y':
        install_dependencies()
    
    # Start the backend in a separate thread
    backend_thread = Thread(target=run_backend)
    backend_thread.daemon = True
    backend_thread.start()
    
    # Give the backend a moment to start
    print("Starting backend server with generative models...")
    time.sleep(2)
    
    # Start the frontend
    print("Starting frontend server...")
    run_frontend()

if __name__ == "__main__":
    main() 