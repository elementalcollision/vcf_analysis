#!/bin/bash
# =============================================================================
# VCF Analysis Agent - Docker Build Script
# =============================================================================
# Comprehensive build script with multi-architecture support, security scanning,
# and flexible build targets for development and production environments.
# =============================================================================

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
IMAGE_NAME="vcf-analysis-agent"
REGISTRY="${REGISTRY:-}"
VERSION="${VERSION:-$(git describe --tags --always --dirty 2>/dev/null || echo 'dev')}"
BUILD_DATE="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
GIT_COMMIT="${GIT_COMMIT:-$(git rev-parse HEAD 2>/dev/null || echo 'unknown')}"

# Build options
TARGET="${TARGET:-runtime}"
PLATFORM="${PLATFORM:-linux/amd64,linux/arm64}"
PUSH="${PUSH:-false}"
SCAN="${SCAN:-true}"
CACHE="${CACHE:-true}"
PARALLEL="${PARALLEL:-true}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Help function
show_help() {
    cat << EOF
VCF Analysis Agent Docker Build Script

Usage: $0 [OPTIONS]

OPTIONS:
    -t, --target TARGET     Build target (runtime|development|builder) [default: runtime]
    -p, --platform PLATFORM Target platforms [default: linux/amd64,linux/arm64]
    -v, --version VERSION   Image version tag [default: git describe]
    -r, --registry REGISTRY Container registry prefix
    --push                  Push images to registry
    --no-cache              Disable build cache
    --no-scan               Skip security scanning
    --no-parallel           Disable parallel builds
    -h, --help              Show this help message

EXAMPLES:
    # Build for local development
    $0 --target development

    # Build and push production image
    $0 --target runtime --push --registry ghcr.io/your-org

    # Build for specific platform
    $0 --platform linux/amd64

    # Build with custom version
    $0 --version v1.2.3

ENVIRONMENT VARIABLES:
    REGISTRY                Container registry prefix
    VERSION                 Image version tag
    GIT_COMMIT             Git commit hash
    DOCKER_BUILDKIT        Enable BuildKit (recommended: 1)
    BUILDX_EXPERIMENTAL    Enable experimental buildx features

EOF
}

# Parse command line arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -t|--target)
                TARGET="$2"
                shift 2
                ;;
            -p|--platform)
                PLATFORM="$2"
                shift 2
                ;;
            -v|--version)
                VERSION="$2"
                shift 2
                ;;
            -r|--registry)
                REGISTRY="$2"
                shift 2
                ;;
            --push)
                PUSH="true"
                shift
                ;;
            --no-cache)
                CACHE="false"
                shift
                ;;
            --no-scan)
                SCAN="false"
                shift
                ;;
            --no-parallel)
                PARALLEL="false"
                shift
                ;;
            -h|--help)
                show_help
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                show_help
                exit 1
                ;;
        esac
    done
}

# Validate environment
validate_environment() {
    log_info "Validating build environment..."
    
    # Check Docker
    if ! command -v docker &> /dev/null; then
        log_error "Docker is not installed or not in PATH"
        exit 1
    fi
    
    # Check Docker Buildx
    if ! docker buildx version &> /dev/null; then
        log_error "Docker Buildx is not available"
        exit 1
    fi
    
    # Enable BuildKit
    export DOCKER_BUILDKIT=1
    export BUILDX_EXPERIMENTAL=1
    
    # Validate target
    case $TARGET in
        runtime|development|builder)
            ;;
        *)
            log_error "Invalid target: $TARGET. Must be one of: runtime, development, builder"
            exit 1
            ;;
    esac
    
    log_success "Environment validation passed"
}

# Setup buildx builder
setup_builder() {
    log_info "Setting up Docker Buildx builder..."
    
    local builder_name="vcf-agent-builder"
    
    # Create builder if it doesn't exist
    if ! docker buildx inspect "$builder_name" &> /dev/null; then
        log_info "Creating new buildx builder: $builder_name"
        docker buildx create \
            --name "$builder_name" \
            --driver docker-container \
            --bootstrap \
            --use
    else
        log_info "Using existing buildx builder: $builder_name"
        docker buildx use "$builder_name"
    fi
    
    # Bootstrap builder
    docker buildx inspect --bootstrap
    
    log_success "Buildx builder ready"
}

# Build image
build_image() {
    log_info "Building VCF Analysis Agent Docker image..."
    log_info "Target: $TARGET"
    log_info "Platform: $PLATFORM"
    log_info "Version: $VERSION"
    
    cd "$PROJECT_ROOT"
    
    # Construct image tag
    local image_tag="$IMAGE_NAME:$VERSION"
    if [[ -n "$REGISTRY" ]]; then
        image_tag="$REGISTRY/$image_tag"
    fi
    
    # Additional tags
    local tags=(
        "--tag" "$image_tag"
        "--tag" "${image_tag}-${TARGET}"
    )
    
    if [[ "$VERSION" != "dev" && "$VERSION" != *"dirty"* ]]; then
        if [[ -n "$REGISTRY" ]]; then
            tags+=("--tag" "$REGISTRY/$IMAGE_NAME:latest")
        else
            tags+=("--tag" "$IMAGE_NAME:latest")
        fi
    fi
    
    # Build arguments
    local build_args=(
        "--build-arg" "VERSION=$VERSION"
        "--build-arg" "BUILD_DATE=$BUILD_DATE"
        "--build-arg" "GIT_COMMIT=$GIT_COMMIT"
        "--build-arg" "BCFTOOLS_VERSION=1.19"
        "--build-arg" "HTSLIB_VERSION=1.19"
    )
    
    # Cache options
    local cache_opts=()
    if [[ "$CACHE" == "true" ]]; then
        cache_opts+=(
            "--cache-from" "type=gha"
            "--cache-to" "type=gha,mode=max"
        )
    fi
    
    # Platform options
    local platform_opts=()
    if [[ "$PARALLEL" == "true" && "$PLATFORM" == *","* ]]; then
        platform_opts+=("--platform" "$PLATFORM")
    else
        # Build for single platform or sequentially
        IFS=',' read -ra PLATFORMS <<< "$PLATFORM"
        for platform in "${PLATFORMS[@]}"; do
            log_info "Building for platform: $platform"
            docker buildx build \
                --target "$TARGET" \
                --platform "$platform" \
                "${tags[@]}" \
                "${build_args[@]}" \
                "${cache_opts[@]}" \
                --load \
                .
        done
        return
    fi
    
    # Output options
    local output_opts=()
    if [[ "$PUSH" == "true" ]]; then
        output_opts+=("--push")
    else
        output_opts+=("--load")
    fi
    
    # Execute build
    log_info "Executing Docker build..."
    docker buildx build \
        --target "$TARGET" \
        "${platform_opts[@]}" \
        "${tags[@]}" \
        "${build_args[@]}" \
        "${cache_opts[@]}" \
        "${output_opts[@]}" \
        .
    
    log_success "Docker build completed successfully"
    
    # Display image information
    if [[ "$PUSH" != "true" ]]; then
        log_info "Built image: $image_tag"
        docker images "$IMAGE_NAME" --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"
    fi
}

# Security scanning
security_scan() {
    if [[ "$SCAN" != "true" ]]; then
        log_info "Skipping security scan (disabled)"
        return
    fi
    
    log_info "Running security scan..."
    
    local image_tag="$IMAGE_NAME:$VERSION"
    if [[ -n "$REGISTRY" ]]; then
        image_tag="$REGISTRY/$image_tag"
    fi
    
    # Try different scanners
    if command -v trivy &> /dev/null; then
        log_info "Scanning with Trivy..."
        trivy image --exit-code 1 --severity HIGH,CRITICAL "$image_tag" || {
            log_warning "Trivy scan found vulnerabilities"
        }
    elif command -v docker &> /dev/null && docker run --rm -v /var/run/docker.sock:/var/run/docker.sock aquasec/trivy:latest --version &> /dev/null; then
        log_info "Scanning with Trivy (Docker)..."
        docker run --rm -v /var/run/docker.sock:/var/run/docker.sock \
            aquasec/trivy:latest image --exit-code 1 --severity HIGH,CRITICAL "$image_tag" || {
            log_warning "Trivy scan found vulnerabilities"
        }
    else
        log_warning "No security scanner available (install trivy for scanning)"
    fi
}

# Cleanup
cleanup() {
    log_info "Cleaning up..."
    
    # Remove dangling images
    docker image prune -f &> /dev/null || true
    
    # Remove build cache if requested
    if [[ "${CLEANUP_CACHE:-false}" == "true" ]]; then
        docker buildx prune -f &> /dev/null || true
    fi
    
    log_success "Cleanup completed"
}

# Main execution
main() {
    log_info "Starting VCF Analysis Agent Docker build..."
    
    parse_args "$@"
    validate_environment
    setup_builder
    build_image
    security_scan
    cleanup
    
    log_success "Build process completed successfully!"
    
    if [[ "$PUSH" != "true" ]]; then
        echo
        log_info "To run the container:"
        echo "  docker run --rm -it $IMAGE_NAME:$VERSION"
        echo
        log_info "To run with docker-compose:"
        echo "  docker-compose up vcf-agent"
    fi
}

# Trap for cleanup on exit
trap cleanup EXIT

# Execute main function
main "$@" 