from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """Settings for the autopvs1 CLI."""

    model_config = SettingsConfigDict(env_file=".env", env_file_encoding="utf-8", case_sensitive=True)

    # --- app-specific settings ---

    # === General settings ===

    DEBUG: bool = False

    # === API settings ===

    #: Base URL to reev
    API_REEV_URL: str = ""


settings = Settings(_env_file=".env", _env_file_encoding="utf-8")  # type: ignore[call-arg]
