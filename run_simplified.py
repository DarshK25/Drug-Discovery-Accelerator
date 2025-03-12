import os
import subprocess
import sys
import time
from threading import Thread

def run_backend():
    """Run the simplified FastAPI backend server."""
    os.chdir('backend')
    subprocess.run([sys.executable, 'simplified_main.py'])

def run_frontend():
    """Run the React frontend development server."""
    os.chdir('frontend')
    subprocess.run(['npm', 'start'])

def main():
    """Run both the frontend and backend servers."""
    print("Starting Drug Discovery Assistant...")
    
    # Start the backend in a separate thread
    backend_thread = Thread(target=run_backend)
    backend_thread.daemon = True
    backend_thread.start()
    
    # Give the backend a moment to start
    print("Starting backend server...")
    time.sleep(2)
    
    # Start the frontend
    print("Starting frontend server...")
    run_frontend()

if __name__ == "__main__":
    main() 