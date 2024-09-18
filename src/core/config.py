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

    #: Whether to enable debug mode
    DEBUG: bool = False
    #: Whether to use cache
    USE_CACHE: bool = True
    #: Path to the cache directory
    CACHE_DIR: str = os.path.join(
        os.path.abspath(os.path.join(__file__, "..", "..", "..")), "cache"
    )

    # === API settings ===

    #: AutoACMG API prefix
    API_V1_STR: str = "/api/v1"

    #: Base URL to reev
    API_REEV_URL: str = ""

    # Keep this setting empty. It is used for the reev-docker-compose.
    #: Base URL to annonars
    AUTO_ACMG_API_ANNONARS_URL: str = ""

    # Keep this setting empty. It is used for the reev-docker-compose.
    #: Base URL to mehari
    AUTO_ACMG_API_MEHARI_URL: str = ""

    # Keep this setting empty. It is used for the reev-docker-compose.
    #: Base URL to dotty
    AUTO_ACMG_API_DOTTY_URL: str = ""

    #: Path to seqrepo data directory
    AUTO_ACMG_SEQREPO_DATA_DIR: str = ""

    #: API key for genebe
    GENEBE_API_KEY: str = ""

    #: Username for genebe
    GENEBE_USERNAME: str = ""

    # ==== Temporary settings ===

    #: Flag to indicate if the duplication is in tandem AND disrupts reading frame AND undergoes NMD
    DUPLICATION_TANDEM: bool = False

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
    api_base_url_annonars: Optional[str] = settings.AUTO_ACMG_API_ANNONARS_URL or None
    #: Base URL to mehari API.
    api_base_url_mehari: Optional[str] = settings.AUTO_ACMG_API_MEHARI_URL or None
    #: Base URL to dotty API.
    api_base_url_dotty: Optional[str] = settings.AUTO_ACMG_API_DOTTY_URL or None

    #: Path to the seqrepo data directory.
    seqrepo_data_dir: Optional[str] = settings.AUTO_ACMG_SEQREPO_DATA_DIR

    @model_validator(mode="before")
    @classmethod
    def _set_base_urls(cls, data: Any) -> Any:
        """If ``api_base_url`` is set, set the other API URLs."""
        if isinstance(data, dict):
            data["api_base_url"] = data.get("api_base_url", settings.API_REEV_URL)
            if (
                not settings.AUTO_ACMG_API_ANNONARS_URL
                and not settings.AUTO_ACMG_API_MEHARI_URL
                and not settings.AUTO_ACMG_API_DOTTY_URL
            ):
                for key in ("annonars", "mehari", "dotty"):
                    data[f"api_base_url_{key}"] = data.get(
                        f"api_base_url_{key}", f"{data['api_base_url']}/{key}"
                    )
        return data
