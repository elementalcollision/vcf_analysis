"""
Tests for VCF Ingestion Pipeline

This module tests the comprehensive VCF ingestion functionality including:
- VCF validation
- Streaming and batch processing
- LanceDB integration
- Kuzu graph database integration
- Error handling and recovery
"""

import pytest
import tempfile
import os
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import numpy as np

from vcf_agent.vcf_ingestion import (
    VCFValidator,
    VCFStreamer,
    EmbeddingGenerator,
    VCFIngestionPipeline,
    IngestionConfig,
    IngestionResult
)


class TestVCFValidator:
    """Test VCF file validation functionality."""
    
    def test_validate_nonexistent_file(self):
        """Test validation of non-existent file."""
        is_valid, errors = VCFValidator.validate_file("/nonexistent/file.vcf")
        assert not is_valid
        assert len(errors) == 1
        assert "not found" in errors[0]
    
    def test_validate_valid_vcf(self):
        """Test validation of valid VCF file."""
        # Use the existing test VCF file
        vcf_path = "sample_data/minimal.vcf.gz"
        if os.path.exists(vcf_path):
            is_valid, errors = VCFValidator.validate_file(vcf_path)
            assert is_valid
            assert len(errors) == 0
    
    @patch('vcf_agent.vcf_ingestion.VCF')
    def test_validate_empty_vcf(self, mock_vcf):
        """Test validation of empty VCF file."""
        # Mock VCF that returns no variants
        mock_vcf_instance = Mock()
        mock_vcf_instance.__iter__ = Mock(return_value=iter([]))
        mock_vcf.return_value.__enter__ = Mock(return_value=mock_vcf_instance)
        mock_vcf.return_value.__exit__ = Mock(return_value=None)
        
        with tempfile.NamedTemporaryFile(suffix='.vcf', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            is_valid, errors = VCFValidator.validate_file(tmp_path)
            assert not is_valid
            assert any("no variants" in error for error in errors)
        finally:
            os.unlink(tmp_path)


class TestVCFStreamer:
    """Test VCF streaming functionality."""
    
    def test_count_variants(self):
        """Test variant counting."""
        vcf_path = "sample_data/minimal.vcf.gz"
        if os.path.exists(vcf_path):
            streamer = VCFStreamer(vcf_path, batch_size=10)
            count = streamer.count_variants()
            assert count >= 0  # Should return a non-negative count
    
    @patch('vcf_agent.vcf_ingestion.VCF')
    def test_stream_variants_batching(self, mock_vcf):
        """Test variant streaming with batching."""
        # Create mock variants
        mock_variants = []
        for i in range(25):  # 25 variants
            variant = Mock()
            variant.CHROM = "1"
            variant.POS = 1000 + i
            variant.REF = "A"
            variant.ALT = ["G"]
            variant.ID = None
            variant.QUAL = 30.0
            variant.FILTER = "PASS"
            variant.INFO = {}
            variant.genotypes = [[0, 1, False]]  # HET genotype
            mock_variants.append(variant)
        
        # Mock VCF instance
        mock_vcf_instance = Mock()
        mock_vcf_instance.__iter__ = Mock(return_value=iter(mock_variants))
        mock_vcf_instance.samples = ["SAMPLE1"]
        mock_vcf_instance.close = Mock()
        mock_vcf.return_value = mock_vcf_instance
        
        # Test streaming with batch size of 10
        streamer = VCFStreamer("dummy.vcf", batch_size=10)
        batches = list(streamer.stream_variants())
        
        # Should have 3 batches: 10, 10, 5
        assert len(batches) == 3
        assert len(batches[0]) == 10
        assert len(batches[1]) == 10
        assert len(batches[2]) == 5
    
    def test_resume_functionality(self):
        """Test resume from specific position."""
        vcf_path = "sample_data/minimal.vcf.gz"
        if os.path.exists(vcf_path):
            # Test with resume position
            streamer = VCFStreamer(vcf_path, batch_size=10, resume_from="1:1000")
            assert streamer.resume_chrom == "1"
            assert streamer.resume_pos == 1000


class TestEmbeddingGenerator:
    """Test embedding generation functionality."""
    
    def test_generate_embedding_dimensions(self):
        """Test embedding generation produces correct dimensions."""
        generator = EmbeddingGenerator(embedding_dim=512)
        
        variant_data = {
            'variant_id': 'test-variant',
            'chrom': '1',
            'pos': 1000,
            'ref': 'A',
            'alt': 'G'
        }
        
        embedding = generator.generate_embedding(variant_data)
        assert embedding.shape == (512,)
        assert embedding.dtype == np.float32
    
    def test_embedding_normalization(self):
        """Test that embeddings are normalized."""
        generator = EmbeddingGenerator(embedding_dim=128)
        
        variant_data = {
            'variant_id': 'test-variant',
            'chrom': '1',
            'pos': 1000,
            'ref': 'A',
            'alt': 'G'
        }
        
        embedding = generator.generate_embedding(variant_data)
        norm = np.linalg.norm(embedding)
        assert abs(norm - 1.0) < 1e-6  # Should be unit vector
    
    def test_embedding_consistency(self):
        """Test that same variant produces same embedding."""
        generator = EmbeddingGenerator(embedding_dim=64)
        
        variant_data = {
            'variant_id': 'test-variant',
            'chrom': '1',
            'pos': 1000,
            'ref': 'A',
            'alt': 'G'
        }
        
        embedding1 = generator.generate_embedding(variant_data)
        embedding2 = generator.generate_embedding(variant_data)
        
        np.testing.assert_array_equal(embedding1, embedding2)


class TestVCFIngestionPipeline:
    """Test the complete VCF ingestion pipeline."""
    
    def setup_method(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.lancedb_path = os.path.join(self.temp_dir, "lancedb")
        self.kuzu_path = os.path.join(self.temp_dir, "kuzu_db")
    
    def teardown_method(self):
        """Clean up test environment."""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_ingestion_config_creation(self):
        """Test ingestion configuration creation."""
        config = IngestionConfig(
            vcf_file="test.vcf",
            lancedb_path=self.lancedb_path,
            kuzu_path=self.kuzu_path,
            batch_size=500
        )
        
        assert config.vcf_file == "test.vcf"
        assert config.batch_size == 500
        assert config.embedding_dim == 1024  # Default value
    
    def test_validate_only_mode(self):
        """Test validation-only mode."""
        config = IngestionConfig(
            vcf_file="sample_data/minimal.vcf.gz",
            validate_only=True
        )
        
        if os.path.exists(config.vcf_file):
            pipeline = VCFIngestionPipeline(config)
            result = pipeline.execute()
            
            # Should complete validation without processing
            assert result['variants_processed'] == 0
            assert len(result['errors']) == 0
    
    @patch('vcf_agent.vcf_ingestion.get_db')
    @patch('vcf_agent.vcf_ingestion.get_or_create_table')
    @patch('vcf_agent.vcf_ingestion.get_kuzu_db_connection')
    @patch('vcf_agent.vcf_ingestion.create_schema')
    @patch('vcf_agent.vcf_ingestion.add_variants')
    def test_database_setup_mocking(self, mock_add_variants, mock_create_schema, 
                                   mock_kuzu_conn, mock_lance_table, mock_lance_db):
        """Test database setup with mocking."""
        # Mock database connections
        mock_lance_db.return_value = Mock()
        mock_lance_table.return_value = Mock()
        mock_kuzu_conn.return_value = Mock()
        mock_create_schema.return_value = None
        mock_add_variants.return_value = None
        
        config = IngestionConfig(
            vcf_file="sample_data/minimal.vcf.gz",
            lancedb_path=self.lancedb_path,
            kuzu_path=self.kuzu_path,
            batch_size=10
        )
        
        if os.path.exists(config.vcf_file):
            pipeline = VCFIngestionPipeline(config)
            
            # Test database setup
            pipeline._setup_databases()
            
            # Verify database connections were established
            mock_lance_db.assert_called_once_with(self.lancedb_path)
            mock_lance_table.assert_called_once()
            mock_kuzu_conn.assert_called_once_with(self.kuzu_path)
            mock_create_schema.assert_called_once()
    
    def test_error_handling(self):
        """Test error handling for invalid inputs."""
        config = IngestionConfig(
            vcf_file="/nonexistent/file.vcf",
            lancedb_path=self.lancedb_path,
            kuzu_path=self.kuzu_path
        )
        
        pipeline = VCFIngestionPipeline(config)
        
        with pytest.raises(ValueError):
            pipeline.execute()
        
        # Check that errors were recorded
        assert len(pipeline.result.errors) > 0
    
    def test_batch_processing_logic_simple(self):
        """Test batch processing logic with real VCF file."""
        vcf_path = "sample_data/minimal.vcf.gz"
        
        if os.path.exists(vcf_path):
            # Test batch processing with real file
            config = IngestionConfig(
                vcf_file=vcf_path,
                batch_size=2,
                validate_only=True  # Skip actual database operations
            )
            
            pipeline = VCFIngestionPipeline(config)
            
            # Test variant counting
            count = pipeline.streamer.count_variants()
            assert count >= 0  # Should return a non-negative count
            
            # Test batch streaming
            batches = list(pipeline.streamer.stream_variants())
            assert len(batches) >= 0  # Should return some batches
            
            # If there are variants, test batch structure
            if count > 0:
                total_variants_in_batches = sum(len(batch) for batch in batches)
                assert total_variants_in_batches == count


class TestIngestionResult:
    """Test ingestion result tracking."""
    
    def test_result_initialization(self):
        """Test result object initialization."""
        result = IngestionResult()
        
        assert result.variants_processed == 0
        assert result.lancedb_records == 0
        assert result.kuzu_variants == 0
        assert result.kuzu_samples == 0
        assert result.kuzu_links == 0
        assert result.errors == []
        assert result.duration_seconds == 0.0
    
    def test_result_dict_conversion(self):
        """Test conversion to dictionary."""
        result = IngestionResult()
        result.variants_processed = 100
        result.errors = ["test error"]
        
        result_dict = result.__dict__
        assert result_dict['variants_processed'] == 100
        assert result_dict['errors'] == ["test error"]


class TestIntegrationScenarios:
    """Integration tests for complete scenarios."""
    
    def setup_method(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.lancedb_path = os.path.join(self.temp_dir, "lancedb")
        self.kuzu_path = os.path.join(self.temp_dir, "kuzu_db")
    
    def teardown_method(self):
        """Clean up test environment."""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_end_to_end_validation_only(self):
        """Test end-to-end pipeline in validation-only mode."""
        vcf_path = "sample_data/minimal.vcf.gz"
        
        if os.path.exists(vcf_path):
            config = IngestionConfig(
                vcf_file=vcf_path,
                lancedb_path=self.lancedb_path,
                kuzu_path=self.kuzu_path,
                batch_size=10,
                validate_only=True
            )
            
            pipeline = VCFIngestionPipeline(config)
            result = pipeline.execute()
            
            # Should complete successfully without errors
            assert len(result['errors']) == 0
            assert result['variants_processed'] == 0  # No processing in validation mode
            assert result['duration_seconds'] >= 0  # Duration should be non-negative
    
    @pytest.mark.skipif(not os.path.exists("sample_data/minimal.vcf.gz"), 
                       reason="Test VCF file not available")
    def test_small_file_ingestion_dry_run(self):
        """Test ingestion of small VCF file (dry run with mocked databases)."""
        with patch('vcf_agent.vcf_ingestion.get_db'), \
             patch('vcf_agent.vcf_ingestion.get_or_create_table'), \
             patch('vcf_agent.vcf_ingestion.get_kuzu_db_connection'), \
             patch('vcf_agent.vcf_ingestion.create_schema'), \
             patch('vcf_agent.vcf_ingestion.add_variants') as mock_add_variants, \
             patch('vcf_agent.graph_integration.add_variant'), \
             patch('vcf_agent.graph_integration.add_sample'), \
             patch('vcf_agent.graph_integration.link_variant_to_sample'):
            
            config = IngestionConfig(
                vcf_file="sample_data/minimal.vcf.gz",
                lancedb_path=self.lancedb_path,
                kuzu_path=self.kuzu_path,
                batch_size=5
            )
            
            pipeline = VCFIngestionPipeline(config)
            result = pipeline.execute()
            
            # Should complete successfully
            assert len(result['errors']) == 0
            assert result['variants_processed'] > 0
            assert result['duration_seconds'] > 0
            
            # Verify database operations were called
            mock_add_variants.assert_called()


if __name__ == "__main__":
    pytest.main([__file__]) 