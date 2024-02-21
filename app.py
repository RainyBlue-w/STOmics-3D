from fastapi import FastAPI
from fastapi.middleware.wsgi import WSGIMiddleware
from dashapp_3dviewer import create_3dviewer
import argparse

# args

parser = argparse.ArgumentParser(description='Run the 3Dviewer app')
parser.add_argument('--port', type=int, default=8000, help='Port to run the app on')
parser.add_argument('--host', type=str, default='0.0.0.0', help='Host to run the app on')
args = parser.parse_args()

# app
app = FastAPI()

@app.get('/')
def root():
    return {'message': 'Hello world from FastAPI!'}

dash_app = create_3dviewer(requests_pathname_prefix='/3dviewer/')
app.mount('/3dviewer/', WSGIMiddleware(dash_app.server))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, port=args.port, host=args.host, log_level="info")
