"""Tests for cofolder.modules.utils.log module."""
import logging


from cofolder.modules.utils.log import setup_root_logger


class TestSetupRootLogger:
    """Tests for setup_root_logger function."""

    def test_setup_with_info_level(self):
        """Test logger setup with INFO level."""
        # Clear any existing handlers first
        logger = logging.getLogger()
        logger.handlers.clear()

        setup_root_logger(logging.INFO)

        # Verify it was configured
        assert logger.level == logging.INFO

    def test_setup_with_debug_level(self):
        """Test logger setup with DEBUG level."""
        # Clear any existing handlers first
        logger = logging.getLogger()
        logger.handlers.clear()

        setup_root_logger(logging.DEBUG)

        # Verify it was configured
        assert logger.level == logging.DEBUG

    def test_setup_with_log_file(self, temp_dir):
        """Test logger setup with log file."""
        # Clear any existing handlers first
        logger = logging.getLogger()
        logger.handlers.clear()

        log_file = temp_dir / "test.log"
        setup_root_logger(logging.INFO, log_file=str(log_file))

        # Log a message
        logger.info("Test message")

        assert log_file.exists()
        content = log_file.read_text()
        assert "Test message" in content

    def test_logger_formatting(self, temp_dir, caplog):
        """Test that logger uses correct format."""
        logger = logging.getLogger()
        logger.handlers.clear()

        log_file = temp_dir / "test.log"
        setup_root_logger(logging.INFO, log_file=str(log_file))

        with caplog.at_level(logging.INFO):
            logger.info("Test message")

        # Check that the message was logged
        assert "Test message" in caplog.text

    def test_reconfigure_preserves_caplog_after_handler_clear(self, temp_dir, caplog):
        """Repeated root logger setup should preserve pytest capture visibility."""
        logger = logging.getLogger()
        logger.handlers.clear()

        setup_root_logger(logging.INFO)

        log_file = temp_dir / "test.log"
        with caplog.at_level(logging.DEBUG):
            logger.handlers.clear()
            setup_root_logger(logging.DEBUG, log_file=str(log_file))
            logger.debug("Debug message after reconfigure")

        assert "Debug message after reconfigure" in caplog.text
        assert log_file.exists()
        assert "Debug message after reconfigure" in log_file.read_text()

    def test_reconfigure_existing_handlers_enables_debug_and_file(self, temp_dir):
        """Repeated setup should upgrade level and add the requested file handler."""
        logger = logging.getLogger()
        logger.handlers.clear()

        setup_root_logger(logging.INFO)

        log_file = temp_dir / "test.log"
        setup_root_logger(logging.DEBUG, log_file=str(log_file))

        logger.debug("Debug message")

        assert logger.level == logging.DEBUG
        assert log_file.exists()
        content = log_file.read_text()
        assert "Debug message" in content
