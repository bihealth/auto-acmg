import os

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
    AUTO_ACMG_API_ANNONARS_URL: str = "http://annonars:8080"

    # Keep this setting empty. It is used for the reev-docker-compose.
    #: Base URL to mehari
    AUTO_ACMG_API_MEHARI_URL: str = "http://mehari:8080"

    # Keep this setting empty. It is used for the reev-docker-compose.
    #: Base URL to dotty
    AUTO_ACMG_API_DOTTY_URL: str = "http://dotty:8080"

    #: Path to seqrepo data directory
    AUTO_ACMG_SEQREPO_DATA_DIR: str = "/home/auto-acmg/seqrepo/master"

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
