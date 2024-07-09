import os
from typing import Any, Optional

from pydantic import BaseModel, ConfigDict, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """Settings for the autopvs1 CLI."""

    model_config = SettingsConfigDict(
        env_file=".env", env_file_encoding="utf-8", case_sensitive=True
    )

    # --- app-specific settings ---

    # === General settings ===

    DEBUG: bool = False

    # === API settings ===

    #: Base URL to reev
    API_REEV_URL: str = ""

    #: Path to seqrepo data directory
    SEQREPO_DATA_DIR: str = ""

    #: API key for genebe
    GENEBE_API_KEY: str = ""

    #: Username for genebe
    GENEBE_USERNAME: str = ""

    #: Path to the root directory
    @property
    def PATH_TO_ROOT(self) -> str:
        return os.path.abspath(os.path.join(__file__, "..", "..", ".."))


settings = Settings(_env_file=".env", _env_file_encoding="utf-8")  # type: ignore[call-arg]


class Config(BaseModel):
    """Configuration for ``Config``."""

    model_config = ConfigDict(frozen=True)

    #: Whether to use dotty for projection (enables legacy transcripts).
    use_dotty: bool = True  # TODO: interpret this in the code!

    #: Base URL to (REEV) API.
    api_base_url: Optional[str] = None
    #: Base URL to annonars API.
    api_base_url_annonars: Optional[str] = None
    #: Base URL to mehari API.
    api_base_url_mehari: Optional[str] = None
    #: Base URL to dotty API.
    api_base_url_dotty: Optional[str] = None

    #: Path to the seqrepo data directory.
    seqrepo_data_dir: Optional[str] = settings.SEQREPO_DATA_DIR

    @model_validator(mode="before")
    @classmethod
    def _set_base_urls(cls, data: Any) -> Any:
        """If ``api_base_url`` is set, set the other API URLs."""
        if isinstance(data, dict):
            data["api_base_url"] = data.get("api_base_url", settings.API_REEV_URL)
            for key in ("annonars", "mehari", "dotty"):
                data[f"api_base_url_{key}"] = data.get(
                    f"api_base_url_{key}", f"{data['api_base_url']}/{key}"
                )
        return data
