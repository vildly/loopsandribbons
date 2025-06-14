from fastapi import FastAPI
from loops_api.routes import router

app = FastAPI(title="LoopsAndRibbons API")

app.include_router(router) 