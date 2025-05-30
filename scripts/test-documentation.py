#!/usr/bin/env python3
"""
Documentation Testing Framework for CLI Enhanced Validation Engine
Comprehensive testing suite for documentation quality, links, performance, and content accuracy.
"""

import subprocess
import pytest
import requests
import time
import json
import re
from pathlib import Path
from typing import List, Dict, Any, Tuple
from urllib.parse import urljoin, urlparse
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class DocumentationTester:
    """Main documentation testing framework."""
    
    def __init__(self, base_dir: Path = None):
        self.base_dir = base_dir or Path(__file__).parent.parent
        self.docs_dir = self.base_dir / "docs" / "docs_src"
        self.site_dir = self.base_dir / "docs_site"
        self.base_url = "http://localhost:8000"
        self.sphinx_build_dir = self.docs_dir / "_build"
        
        # Performance thresholds
        self.performance_thresholds = {
            'build_time_mkdocs': 10.0,  # seconds
            'build_time_sphinx': 15.0,  # seconds
            'page_load_time': 3.0,      # seconds
            'link_response_time': 5.0,   # seconds
            'total_build_time': 20.0     # seconds
        }
        
        # Link validation settings
        self.external_link_timeout = 10
        self.max_concurrent_links = 10
        
    def run_all_tests(self) -> Dict[str, Any]:
        """Run all documentation tests and return comprehensive results."""
        logger.info("🧪 Starting comprehensive documentation testing...")
        
        results = {
            'start_time': time.time(),
            'tests': {},
            'performance_metrics': {},
            'summary': {}
        }
        
        try:
            # Phase 1: Build Tests
            results['tests']['mkdocs_build'] = self.test_mkdocs_build()
            results['tests']['sphinx_build'] = self.test_sphinx_build()
            
            # Phase 2: Content Tests
            results['tests']['content_validation'] = self.test_content_validation()
            results['tests']['link_validation'] = self.test_link_validation()
            
            # Phase 3: Performance Tests
            results['tests']['performance'] = self.test_performance()
            
            # Phase 4: Quality Tests
            results['tests']['accessibility'] = self.test_accessibility()
            results['tests']['seo'] = self.test_seo()
            
        except Exception as e:
            logger.error(f"Testing failed: {e}")
            results['error'] = str(e)
        
        finally:
            results['end_time'] = time.time()
            results['total_duration'] = results['end_time'] - results['start_time']
            results['summary'] = self._generate_summary(results)
            
        return results
    
    def test_mkdocs_build(self) -> Dict[str, Any]:
        """Test MkDocs builds successfully without errors."""
        logger.info("🏗️ Testing MkDocs build...")
        
        start_time = time.time()
        
        try:
            # Clean previous build
            if self.site_dir.exists():
                import shutil
                shutil.rmtree(self.site_dir)
            
            # Run MkDocs build
            result = subprocess.run(
                ['mkdocs', 'build', '--strict'],
                cwd=self.base_dir,
                capture_output=True,
                text=True,
                timeout=self.performance_thresholds['build_time_mkdocs']
            )
            
            build_time = time.time() - start_time
            
            return {
                'success': result.returncode == 0,
                'build_time': build_time,
                'stdout': result.stdout,
                'stderr': result.stderr,
                'files_generated': self._count_generated_files(self.site_dir),
                'performance_ok': build_time < self.performance_thresholds['build_time_mkdocs']
            }
            
        except subprocess.TimeoutExpired:
            return {
                'success': False,
                'error': f"Build timeout after {self.performance_thresholds['build_time_mkdocs']}s",
                'build_time': self.performance_thresholds['build_time_mkdocs']
            }
        except Exception as e:
            return {
                'success': False,
                'error': str(e),
                'build_time': time.time() - start_time
            }
    
    def test_sphinx_build(self) -> Dict[str, Any]:
        """Test Sphinx API documentation generates correctly."""
        logger.info("📚 Testing Sphinx build...")
        
        start_time = time.time()
        
        try:
            # Clean previous build
            if self.sphinx_build_dir.exists():
                import shutil
                shutil.rmtree(self.sphinx_build_dir)
            
            # Run Sphinx build
            result = subprocess.run(
                ['sphinx-build', '-b', 'html', '.', '_build', '-W'],
                cwd=self.docs_dir,
                capture_output=True,
                text=True,
                timeout=self.performance_thresholds['build_time_sphinx']
            )
            
            build_time = time.time() - start_time
            
            return {
                'success': result.returncode == 0,
                'build_time': build_time,
                'stdout': result.stdout,
                'stderr': result.stderr,
                'api_files_generated': self._count_generated_files(self.sphinx_build_dir),
                'performance_ok': build_time < self.performance_thresholds['build_time_sphinx']
            }
            
        except subprocess.TimeoutExpired:
            return {
                'success': False,
                'error': f"Sphinx build timeout after {self.performance_thresholds['build_time_sphinx']}s",
                'build_time': self.performance_thresholds['build_time_sphinx']
            }
        except Exception as e:
            return {
                'success': False,
                'error': str(e),
                'build_time': time.time() - start_time
            }
    
    def test_content_validation(self) -> Dict[str, Any]:
        """Validate content structure and completeness."""
        logger.info("📝 Testing content validation...")
        
        issues = []
        checks = {
            'required_files': 0,
            'broken_links': 0,
            'missing_images': 0,
            'malformed_markdown': 0
        }
        
        # Check required files exist
        required_files = [
            'index.md',
            'user_guide/index.md',
            'user_guide/quick_start.md',
            'user_guide/configuration.md',
            'user_guide/cli_reference.md',
            'user_guide/examples.md',
            'user_guide/troubleshooting.md',
            'integration/index.md',
            'integration/github_actions.md',
            'api/index.md'
        ]
        
        for file_path in required_files:
            full_path = self.docs_dir / file_path
            if not full_path.exists():
                issues.append(f"Missing required file: {file_path}")
                checks['required_files'] += 1
        
        # Validate markdown files
        for md_file in self.docs_dir.rglob("*.md"):
            try:
                content = md_file.read_text(encoding='utf-8')
                
                # Check for broken internal links
                internal_links = re.findall(r'\[([^\]]+)\]\(([^)]+)\)', content)
                for link_text, link_url in internal_links:
                    if link_url.startswith('../') or link_url.startswith('./'):
                        target_path = (md_file.parent / link_url).resolve()
                        if not target_path.exists() and not target_path.with_suffix('.md').exists():
                            issues.append(f"Broken internal link in {md_file.relative_to(self.docs_dir)}: {link_url}")
                            checks['broken_links'] += 1
                
                # Check for missing images
                image_links = re.findall(r'!\[([^\]]*)\]\(([^)]+)\)', content)
                for alt_text, img_url in image_links:
                    if not img_url.startswith('http'):
                        img_path = (md_file.parent / img_url).resolve()
                        if not img_path.exists():
                            issues.append(f"Missing image in {md_file.relative_to(self.docs_dir)}: {img_url}")
                            checks['missing_images'] += 1
                
            except Exception as e:
                issues.append(f"Error reading {md_file.relative_to(self.docs_dir)}: {e}")
                checks['malformed_markdown'] += 1
        
        return {
            'success': len(issues) == 0,
            'issues': issues,
            'checks': checks,
            'total_issues': len(issues)
        }
    
    def test_link_validation(self) -> Dict[str, Any]:
        """Validate all internal and external links."""
        logger.info("🔗 Testing link validation...")
        
        if not self.site_dir.exists():
            return {
                'success': False,
                'error': 'Site directory does not exist. Run build first.'
            }
        
        all_links = self._extract_all_links()
        internal_links = [link for link in all_links if self._is_internal_link(link)]
        external_links = [link for link in all_links if not self._is_internal_link(link)]
        
        # Test internal links
        internal_results = self._test_internal_links(internal_links)
        
        # Test external links (with timeout and concurrency)
        external_results = self._test_external_links(external_links)
        
        total_broken = internal_results['broken'] + external_results['broken']
        
        return {
            'success': total_broken == 0,
            'internal_links': internal_results,
            'external_links': external_results,
            'total_links': len(all_links),
            'total_broken': total_broken,
            'broken_percentage': (total_broken / len(all_links) * 100) if all_links else 0
        }
    
    def test_performance(self) -> Dict[str, Any]:
        """Test page load performance and build times."""
        logger.info("⚡ Testing performance...")
        
        # Start a local server for testing
        server_process = None
        try:
            server_process = subprocess.Popen(
                ['mkdocs', 'serve', '--dev-addr', '127.0.0.1:8000'],
                cwd=self.base_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            
            # Wait for server to start
            time.sleep(3)
            
            # Test key pages
            test_pages = [
                '/',
                '/user_guide/quick_start/',
                '/user_guide/cli_reference/',
                '/integration/github_actions/',
                '/api/'
            ]
            
            page_results = {}
            for page in test_pages:
                page_results[page] = self._test_page_performance(page)
            
            return {
                'success': all(result['success'] for result in page_results.values()),
                'pages': page_results,
                'average_load_time': sum(r.get('load_time', 0) for r in page_results.values()) / len(page_results),
                'performance_threshold': self.performance_thresholds['page_load_time']
            }
            
        except Exception as e:
            return {
                'success': False,
                'error': str(e)
            }
        finally:
            if server_process:
                server_process.terminate()
                server_process.wait()
    
    def test_accessibility(self) -> Dict[str, Any]:
        """Test accessibility compliance."""
        logger.info("♿ Testing accessibility...")
        
        issues = []
        checks = {
            'missing_alt_text': 0,
            'missing_headings': 0,
            'poor_contrast': 0,
            'missing_landmarks': 0
        }
        
        # Check HTML files for accessibility issues
        for html_file in self.site_dir.rglob("*.html"):
            try:
                content = html_file.read_text(encoding='utf-8')
                
                # Check for images without alt text
                img_without_alt = re.findall(r'<img(?![^>]*alt=)[^>]*>', content)
                if img_without_alt:
                    issues.append(f"Images without alt text in {html_file.name}: {len(img_without_alt)}")
                    checks['missing_alt_text'] += len(img_without_alt)
                
                # Check for proper heading structure
                headings = re.findall(r'<h([1-6])[^>]*>', content)
                if headings:
                    # Check if headings follow proper order
                    heading_levels = [int(h) for h in headings]
                    for i in range(1, len(heading_levels)):
                        if heading_levels[i] - heading_levels[i-1] > 1:
                            issues.append(f"Heading level skip in {html_file.name}")
                            checks['missing_headings'] += 1
                            break
                
                # Check for main landmarks
                if '<main' not in content and 'role="main"' not in content:
                    issues.append(f"Missing main landmark in {html_file.name}")
                    checks['missing_landmarks'] += 1
                
            except Exception as e:
                issues.append(f"Error checking accessibility in {html_file.name}: {e}")
        
        return {
            'success': len(issues) == 0,
            'issues': issues,
            'checks': checks,
            'total_issues': len(issues)
        }
    
    def test_seo(self) -> Dict[str, Any]:
        """Test SEO optimization."""
        logger.info("🔍 Testing SEO optimization...")
        
        issues = []
        checks = {
            'missing_title': 0,
            'missing_description': 0,
            'missing_og_tags': 0,
            'duplicate_titles': 0
        }
        
        titles_found = []
        
        for html_file in self.site_dir.rglob("*.html"):
            try:
                content = html_file.read_text(encoding='utf-8')
                
                # Check for title tag
                title_match = re.search(r'<title[^>]*>([^<]+)</title>', content)
                if not title_match:
                    issues.append(f"Missing title tag in {html_file.name}")
                    checks['missing_title'] += 1
                else:
                    title = title_match.group(1).strip()
                    if title in titles_found:
                        issues.append(f"Duplicate title in {html_file.name}: {title}")
                        checks['duplicate_titles'] += 1
                    titles_found.append(title)
                
                # Check for meta description
                if 'name="description"' not in content:
                    issues.append(f"Missing meta description in {html_file.name}")
                    checks['missing_description'] += 1
                
                # Check for Open Graph tags
                if 'property="og:' not in content:
                    issues.append(f"Missing Open Graph tags in {html_file.name}")
                    checks['missing_og_tags'] += 1
                
            except Exception as e:
                issues.append(f"Error checking SEO in {html_file.name}: {e}")
        
        return {
            'success': len(issues) == 0,
            'issues': issues,
            'checks': checks,
            'total_issues': len(issues)
        }
    
    def _count_generated_files(self, directory: Path) -> int:
        """Count generated files in a directory."""
        if not directory.exists():
            return 0
        return len(list(directory.rglob("*")))
    
    def _extract_all_links(self) -> List[str]:
        """Extract all links from generated HTML files."""
        links = set()
        
        for html_file in self.site_dir.rglob("*.html"):
            try:
                content = html_file.read_text(encoding='utf-8')
                # Extract href attributes
                href_matches = re.findall(r'href=["\']([^"\']+)["\']', content)
                links.update(href_matches)
            except Exception as e:
                logger.warning(f"Error extracting links from {html_file}: {e}")
        
        return list(links)
    
    def _is_internal_link(self, link: str) -> bool:
        """Check if a link is internal."""
        return (
            link.startswith('/') or 
            link.startswith('./') or 
            link.startswith('../') or
            link.startswith('#') or
            not link.startswith('http')
        )
    
    def _test_internal_links(self, links: List[str]) -> Dict[str, Any]:
        """Test internal links."""
        broken_links = []
        
        for link in links:
            if link.startswith('#'):
                continue  # Skip anchors for now
            
            # Convert to file path
            if link.startswith('/'):
                file_path = self.site_dir / link.lstrip('/')
            else:
                file_path = self.site_dir / link
            
            # Check if it's a directory (should have index.html)
            if file_path.is_dir():
                file_path = file_path / 'index.html'
            elif not file_path.suffix:
                file_path = file_path.with_suffix('.html')
            
            if not file_path.exists():
                broken_links.append(link)
        
        return {
            'total': len(links),
            'broken': len(broken_links),
            'broken_links': broken_links,
            'success_rate': ((len(links) - len(broken_links)) / len(links) * 100) if links else 100
        }
    
    def _test_external_links(self, links: List[str]) -> Dict[str, Any]:
        """Test external links with concurrency and timeout."""
        broken_links = []
        
        def test_link(link: str) -> Tuple[str, bool]:
            try:
                response = requests.head(
                    link, 
                    timeout=self.external_link_timeout,
                    allow_redirects=True
                )
                return link, response.status_code < 400
            except Exception:
                return link, False
        
        # Test external links concurrently
        with ThreadPoolExecutor(max_workers=self.max_concurrent_links) as executor:
            future_to_link = {
                executor.submit(test_link, link): link 
                for link in links if link.startswith('http')
            }
            
            for future in as_completed(future_to_link):
                link, is_valid = future.result()
                if not is_valid:
                    broken_links.append(link)
        
        return {
            'total': len([l for l in links if l.startswith('http')]),
            'broken': len(broken_links),
            'broken_links': broken_links,
            'success_rate': ((len(links) - len(broken_links)) / len(links) * 100) if links else 100
        }
    
    def _test_page_performance(self, page: str) -> Dict[str, Any]:
        """Test individual page performance."""
        url = urljoin(self.base_url, page)
        
        try:
            start_time = time.time()
            response = requests.get(url, timeout=self.performance_thresholds['page_load_time'])
            load_time = time.time() - start_time
            
            return {
                'success': response.status_code == 200 and load_time < self.performance_thresholds['page_load_time'],
                'status_code': response.status_code,
                'load_time': load_time,
                'content_length': len(response.content),
                'performance_ok': load_time < self.performance_thresholds['page_load_time']
            }
            
        except Exception as e:
            return {
                'success': False,
                'error': str(e),
                'load_time': self.performance_thresholds['page_load_time']
            }
    
    def _generate_summary(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Generate test summary."""
        total_tests = len(results['tests'])
        passed_tests = sum(1 for test in results['tests'].values() if test.get('success', False))
        
        return {
            'total_tests': total_tests,
            'passed_tests': passed_tests,
            'failed_tests': total_tests - passed_tests,
            'success_rate': (passed_tests / total_tests * 100) if total_tests else 0,
            'total_duration': results.get('total_duration', 0),
            'overall_success': passed_tests == total_tests
        }


def main():
    """Main test runner."""
    tester = DocumentationTester()
    results = tester.run_all_tests()
    
    # Print results
    print(json.dumps(results, indent=2, default=str))
    
    # Exit with appropriate code
    exit_code = 0 if results['summary']['overall_success'] else 1
    return exit_code


if __name__ == "__main__":
    exit(main()) 