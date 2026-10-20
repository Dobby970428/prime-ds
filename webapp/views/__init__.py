"""블루프린트 등록."""

from __future__ import annotations

from flask import Flask


def register_blueprints(app: Flask) -> None:
    from .api import bp as api_bp
    from .degradation import bp as degradation_bp
    from .excipient import bp as excipient_bp
    from .main import bp as main_bp
    from .process import bp as process_bp
    from .reports import bp as reports_bp
    from .search import bp as search_bp
    from .single import bp as single_bp

    app.register_blueprint(main_bp)
    app.register_blueprint(single_bp, url_prefix="/single")
    app.register_blueprint(process_bp, url_prefix="/process")
    app.register_blueprint(degradation_bp, url_prefix="/degradation")
    app.register_blueprint(excipient_bp, url_prefix="/excipient")
    app.register_blueprint(search_bp, url_prefix="/search")
    app.register_blueprint(reports_bp, url_prefix="/reports")
    app.register_blueprint(api_bp, url_prefix="/api")

    from .context import register_context

    register_context(app)
