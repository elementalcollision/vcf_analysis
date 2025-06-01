/**
 * Performance Optimization JavaScript for CLI Enhanced Validation Engine Documentation
 * Features: Performance monitoring, lazy loading, prefetching, and UX enhancements
 */

(function() {
    'use strict';

    // Performance monitoring
    const Performance = {
        metrics: {},
        
        init() {
            this.measurePageLoad();
            this.setupPerformanceObserver();
            this.monitorResourceLoading();
            this.trackUserInteractions();
        },

        measurePageLoad() {
            if (window.performance && window.performance.timing) {
                const timing = window.performance.timing;
                const loadTime = timing.loadEventEnd - timing.navigationStart;
                this.metrics.pageLoad = loadTime;
                
                // Show performance indicator in debug mode
                if (window.location.search.includes('debug=performance')) {
                    this.showPerformanceIndicator(loadTime);
                }
            }
        },

        setupPerformanceObserver() {
            if ('PerformanceObserver' in window) {
                // Observe Largest Contentful Paint
                const lcpObserver = new PerformanceObserver((list) => {
                    const entries = list.getEntries();
                    const lastEntry = entries[entries.length - 1];
                    this.metrics.lcp = lastEntry.startTime;
                });
                lcpObserver.observe({ entryTypes: ['largest-contentful-paint'] });

                // Observe First Input Delay
                const fidObserver = new PerformanceObserver((list) => {
                    const entries = list.getEntries();
                    entries.forEach((entry) => {
                        this.metrics.fid = entry.processingStart - entry.startTime;
                    });
                });
                fidObserver.observe({ entryTypes: ['first-input'] });

                // Observe Cumulative Layout Shift
                const clsObserver = new PerformanceObserver((list) => {
                    const entries = list.getEntries();
                    let cls = 0;
                    entries.forEach((entry) => {
                        if (!entry.hadRecentInput) {
                            cls += entry.value;
                        }
                    });
                    this.metrics.cls = cls;
                });
                clsObserver.observe({ entryTypes: ['layout-shift'] });
            }
        },

        monitorResourceLoading() {
            if ('PerformanceObserver' in window) {
                const resourceObserver = new PerformanceObserver((list) => {
                    const entries = list.getEntries();
                    entries.forEach((entry) => {
                        if (entry.duration > 1000) { // Resources taking >1s
                            console.warn(`Slow resource: ${entry.name} took ${entry.duration}ms`);
                        }
                    });
                });
                resourceObserver.observe({ entryTypes: ['resource'] });
            }
        },

        trackUserInteractions() {
            let interactionCount = 0;
            const trackInteraction = () => {
                interactionCount++;
                if (interactionCount === 1) {
                    this.metrics.timeToFirstInteraction = performance.now();
                }
            };

            ['click', 'keydown', 'scroll'].forEach(event => {
                document.addEventListener(event, trackInteraction, { once: true, passive: true });
            });
        },

        showPerformanceIndicator(loadTime) {
            const indicator = document.createElement('div');
            indicator.className = 'performance-indicator';
            indicator.textContent = `Load: ${loadTime}ms`;
            document.body.appendChild(indicator);
            
            // Update with other metrics as they become available
            setTimeout(() => {
                let content = `Load: ${loadTime}ms`;
                if (this.metrics.lcp) content += ` | LCP: ${Math.round(this.metrics.lcp)}ms`;
                if (this.metrics.fid) content += ` | FID: ${Math.round(this.metrics.fid)}ms`;
                if (this.metrics.cls) content += ` | CLS: ${this.metrics.cls.toFixed(3)}`;
                indicator.textContent = content;
            }, 3000);
        }
    };

    // Enhanced lazy loading for images
    const LazyLoading = {
        init() {
            if ('IntersectionObserver' in window) {
                this.setupImageLazyLoading();
                this.setupContentLazyLoading();
            }
        },

        setupImageLazyLoading() {
            const imageObserver = new IntersectionObserver((entries) => {
                entries.forEach(entry => {
                    if (entry.isIntersecting) {
                        const img = entry.target;
                        if (img.dataset.src) {
                            img.src = img.dataset.src;
                            img.removeAttribute('data-src');
                            imageObserver.unobserve(img);
                        }
                    }
                });
            }, { rootMargin: '50px' });

            document.querySelectorAll('img[data-src]').forEach(img => {
                imageObserver.observe(img);
            });
        },

        setupContentLazyLoading() {
            const contentObserver = new IntersectionObserver((entries) => {
                entries.forEach(entry => {
                    if (entry.isIntersecting) {
                        const element = entry.target;
                        element.classList.add('lazy-loaded');
                        
                        // Load syntax highlighting for code blocks
                        if (element.classList.contains('highlight')) {
                            this.loadSyntaxHighlighting(element);
                        }
                        
                        contentObserver.unobserve(element);
                    }
                });
            }, { rootMargin: '100px' });

            document.querySelectorAll('.lazy-content').forEach(element => {
                contentObserver.observe(element);
            });
        },

        loadSyntaxHighlighting(element) {
            // Enhance syntax highlighting for performance
            const codeBlocks = element.querySelectorAll('code');
            codeBlocks.forEach(block => {
                if (!block.classList.contains('hljs-processed')) {
                    block.classList.add('hljs-processed');
                }
            });
        }
    };

    // Link prefetching for instant navigation
    const Prefetching = {
        prefetchedUrls: new Set(),

        init() {
            this.setupLinkPrefetching();
            this.setupPredictivePrefetching();
        },

        setupLinkPrefetching() {
            if ('IntersectionObserver' in window) {
                const linkObserver = new IntersectionObserver((entries) => {
                    entries.forEach(entry => {
                        if (entry.isIntersecting) {
                            const link = entry.target;
                            this.prefetchLink(link.href);
                            linkObserver.unobserve(link);
                        }
                    });
                }, { rootMargin: '200px' });

                document.querySelectorAll('a[href^="/"], a[href^="./"], a[href^="../"]').forEach(link => {
                    linkObserver.observe(link);
                });
            }
        },

        setupPredictivePrefetching() {
            let hoverTimer;
            
            document.addEventListener('mouseover', (e) => {
                if (e.target.tagName === 'A' && e.target.href) {
                    hoverTimer = setTimeout(() => {
                        this.prefetchLink(e.target.href);
                    }, 100);
                }
            });

            document.addEventListener('mouseout', () => {
                if (hoverTimer) {
                    clearTimeout(hoverTimer);
                }
            });
        },

        prefetchLink(url) {
            if (this.prefetchedUrls.has(url)) return;
            
            const link = document.createElement('link');
            link.rel = 'prefetch';
            link.href = url;
            document.head.appendChild(link);
            
            this.prefetchedUrls.add(url);
        }
    };

    // Progressive enhancement utilities
    const ProgressiveEnhancement = {
        init() {
            this.enhanceSearch();
            this.enhanceNavigation();
            this.enhanceCodeBlocks();
            this.enhanceAccessibility();
        },

        enhanceSearch() {
            const searchInput = document.querySelector('.md-search__input');
            if (searchInput) {
                // Debounce search input
                let searchTimeout;
                searchInput.addEventListener('input', (e) => {
                    clearTimeout(searchTimeout);
                    searchTimeout = setTimeout(() => {
                        this.performSearch(e.target.value);
                    }, 300);
                });
            }
        },

        performSearch(query) {
            // Enhanced search functionality
            if (query.length > 2) {
                // Implement enhanced search logic here
                console.log('Searching for:', query);
            }
        },

        enhanceNavigation() {
            // Add smooth scrolling for anchor links
            document.querySelectorAll('a[href^="#"]').forEach(anchor => {
                anchor.addEventListener('click', (e) => {
                    e.preventDefault();
                    const target = document.querySelector(anchor.getAttribute('href'));
                    if (target) {
                        target.scrollIntoView({
                            behavior: 'smooth',
                            block: 'start'
                        });
                    }
                });
            });

            // Add active navigation highlighting
            this.updateActiveNavigation();
            window.addEventListener('scroll', () => {
                this.throttle(this.updateActiveNavigation.bind(this), 100)();
            });
        },

        updateActiveNavigation() {
            const sections = document.querySelectorAll('[id]');
            const navLinks = document.querySelectorAll('.md-nav__link');
            
            let current = '';
            sections.forEach(section => {
                const rect = section.getBoundingClientRect();
                if (rect.top <= 100) {
                    current = section.id;
                }
            });

            navLinks.forEach(link => {
                link.classList.remove('md-nav__link--active');
                if (link.getAttribute('href') === `#${current}`) {
                    link.classList.add('md-nav__link--active');
                }
            });
        },

        enhanceCodeBlocks() {
            document.querySelectorAll('.highlight').forEach(block => {
                // Add copy button if not already present
                if (!block.querySelector('.md-clipboard')) {
                    this.addCopyButton(block);
                }

                // Add line numbers if requested
                if (block.dataset.lineNumbers) {
                    this.addLineNumbers(block);
                }
            });
        },

        addCopyButton(codeBlock) {
            const button = document.createElement('button');
            button.className = 'md-clipboard md-icon';
            button.title = 'Copy to clipboard';
            button.innerHTML = '📋';
            
            button.addEventListener('click', () => {
                const code = codeBlock.querySelector('code');
                if (code) {
                    navigator.clipboard.writeText(code.textContent).then(() => {
                        button.innerHTML = '✅';
                        setTimeout(() => {
                            button.innerHTML = '📋';
                        }, 2000);
                    });
                }
            });

            codeBlock.style.position = 'relative';
            codeBlock.appendChild(button);
        },

        addLineNumbers(codeBlock) {
            const code = codeBlock.querySelector('code');
            if (code) {
                const lines = code.textContent.split('\n').length;
                const lineNumbers = document.createElement('div');
                lineNumbers.className = 'line-numbers';
                
                for (let i = 1; i <= lines; i++) {
                    const line = document.createElement('span');
                    line.textContent = i;
                    lineNumbers.appendChild(line);
                }
                
                codeBlock.insertBefore(lineNumbers, code);
            }
        },

        enhanceAccessibility() {
            // Add skip to content link
            const skipLink = document.createElement('a');
            skipLink.href = '#main-content';
            skipLink.textContent = 'Skip to main content';
            skipLink.className = 'skip-link';
            document.body.insertBefore(skipLink, document.body.firstChild);

            // Enhance focus management
            document.addEventListener('keydown', (e) => {
                if (e.key === 'Tab') {
                    document.body.classList.add('keyboard-navigation');
                }
            });

            document.addEventListener('mousedown', () => {
                document.body.classList.remove('keyboard-navigation');
            });
        },

        throttle(func, limit) {
            let inThrottle;
            return function() {
                const args = arguments;
                const context = this;
                if (!inThrottle) {
                    func.apply(context, args);
                    inThrottle = true;
                    setTimeout(() => inThrottle = false, limit);
                }
            };
        }
    };

    // Error tracking and reporting
    const ErrorTracking = {
        init() {
            window.addEventListener('error', (e) => {
                this.logError('JavaScript Error', e.error);
            });

            window.addEventListener('unhandledrejection', (e) => {
                this.logError('Unhandled Promise Rejection', e.reason);
            });
        },

        logError(type, error) {
            console.error(`${type}:`, error);
            
            // In production, send to error tracking service
            if (window.location.hostname !== 'localhost') {
                // Example: Send to error tracking service
                // this.sendToErrorService(type, error);
            }
        }
    };

    // Initialize all performance features
    /**
    * Initializes various performance optimization features when the DOM is ready.
    * @example
    * init()
    * // Logs "📊 Performance optimization features initialized"
    * @param {void} None - This function does not take any parameters.
    * @returns {void} This function does not return anything.
    * @description
    *   - Sets up event listener to delay initialization until the DOM is fully loaded.
    *   - Initializes features like Performance, LazyLoading, Prefetching, ProgressiveEnhancement, and ErrorTracking.
    *   - Adds a loading progress indicator to the body of the document.
    *   - Ensures progress bar activation during page navigation.
    */
    function init() {
        // Wait for DOM to be ready
        if (document.readyState === 'loading') {
            document.addEventListener('DOMContentLoaded', init);
            return;
        }

        Performance.init();
        LazyLoading.init();
        Prefetching.init();
        ProgressiveEnhancement.init();
        ErrorTracking.init();

        // Add loading progress indicator
        const progressBar = document.createElement('div');
        progressBar.className = 'loading-progress';
        document.body.appendChild(progressBar);

        // Show progress during navigation
        window.addEventListener('beforeunload', () => {
            progressBar.classList.add('active');
        });

        console.log('📊 Performance optimization features initialized');
    }

    // Initialize when script loads
    init();

    // Expose performance API for debugging
    window.DocumentationPerformance = {
        metrics: Performance.metrics,
        prefetched: Prefetching.prefetchedUrls
    };
})(); 