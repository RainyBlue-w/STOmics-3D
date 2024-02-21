from fastapi import FastAPI
from fastapi.middleware.wsgi import WSGIMiddleware
from dashapp_3dviewer import create_3dviewer

app = FastAPI()

@app.get('/')
def root():
    return {'message': 'Hello world from FastAPI!'}

dash_app = create_3dviewer(requests_pathname_prefix='/dash/')
app.mount('/dash', WSGIMiddleware(dash_app.server))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, port=8051)