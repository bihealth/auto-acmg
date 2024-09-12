"""Entry point for the AutoACMG API."""

import pathlib

from fastapi import FastAPI
from fastapi.responses import FileResponse

from src.api.internal.api import router as internal_router
from src.core.config import settings

app = FastAPI(
    title="AutoACMG API",
    description="API for AutoACMG",
    version="1.0.0",
    docs_url=f"{settings.API_V1_STR}/docs",
    openapi_url=f"{settings.API_V1_STR}/openapi.json",
    debug=settings.DEBUG,
)


@app.get("/favicon.ico", include_in_schema=False)
async def favicon():
    """Serve favicon"""
    return FileResponse(pathlib.Path(__file__).parent / "assets/favicon.ico")


app.include_router(internal_router, prefix=f"{settings.API_V1_STR}/internal")
