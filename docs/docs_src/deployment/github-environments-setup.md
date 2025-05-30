# GitHub Environments Setup Guide

**Last Updated**: May 29, 2025  
**Required for**: Production deployment workflow (`.github/workflows/production-deploy.yml`)

## Overview

This guide explains how to set up GitHub environments required for the VCF Analysis Agent production deployment workflow. Environments provide deployment protection, approval processes, and environment-specific secrets management.

## Prerequisites

- Repository admin access
- GitHub Pro, Team, or Enterprise plan (for private repositories)
- Understanding of GitHub Actions and deployment workflows

## Environment Configuration

### 1. Staging Environment

**Purpose**: Pre-production testing and validation  
**URL**: `https://vcf-agent-staging.example.com`  
**Protection Level**: Minimal (for rapid iteration)

#### Setup Steps:
1. Navigate to repository **Settings** > **Environments**
2. Click **New environment**
3. Name: `staging`
4. Configure settings:
   ```yaml
   Environment name: staging
   Deployment branches: Any branch
   Required reviewers: None (optional: add 1 reviewer)
   Wait timer: 0 minutes
   Environment URL: https://vcf-agent-staging.example.com
   ```

#### Required Secrets:
- `SLACK_WEBHOOK_URL`: Notification webhook for deployment status
- `STAGING_DATABASE_URL`: Database connection for staging
- `STAGING_API_KEYS`: Third-party service API keys

#### Environment Variables:
- `ENVIRONMENT=staging`
- `LOG_LEVEL=DEBUG`
- `OTEL_SAMPLING_RATE=1.0`

### 2. Production Environment

**Purpose**: Live production deployment  
**URL**: `https://vcf-agent.example.com`  
**Protection Level**: Maximum (multi-reviewer approval required)

#### Setup Steps:
1. Navigate to repository **Settings** > **Environments**
2. Click **New environment**
3. Name: `production`
4. Configure settings:
   ```yaml
   Environment name: production
   Deployment branches: 
     - main
     - tags matching pattern: v*
   Required reviewers: 
     - @production-team
     - @security-team
   Wait timer: 5 minutes (cooling-off period)
   Environment URL: https://vcf-agent.example.com
   ```

#### Protection Rules:
- **Required reviewers**: Minimum 2 reviewers from production team
- **Deployment branches**: Only `main` branch and version tags (`v*`)
- **Wait timer**: 5-minute mandatory delay for critical assessments
- **Prevent self-review**: Deployer cannot approve their own deployment

#### Required Secrets:
- `AWS_ACCESS_KEY_ID`: Production AWS access key
- `AWS_SECRET_ACCESS_KEY`: Production AWS secret key
- `SLACK_WEBHOOK_URL`: Production notification webhook
- `PRODUCTION_DATABASE_URL`: Production database connection
- `PRODUCTION_API_KEYS`: Production API keys and certificates
- `NEWRELIC_LICENSE_KEY`: Production monitoring (optional)
- `DATADOG_API_KEY`: Production observability (optional)

#### Environment Variables:
- `ENVIRONMENT=production`
- `LOG_LEVEL=INFO`
- `OTEL_SAMPLING_RATE=0.1`
- `BACKUP_ENABLED=true`
- `HIGH_AVAILABILITY=true`

### 3. Production Rollback Environment

**Purpose**: Emergency rollback capability  
**URL**: `https://vcf-agent.example.com`  
**Protection Level**: Emergency access (limited reviewers)

#### Setup Steps:
1. Navigate to repository **Settings** > **Environments**
2. Click **New environment** 
3. Name: `production-rollback`
4. Configure settings:
   ```yaml
   Environment name: production-rollback
   Deployment branches: Any branch (emergency access)
   Required reviewers: 
     - @incident-response-team
     - @sre-team
   Wait timer: 0 minutes (emergency response)
   Environment URL: https://vcf-agent.example.com
   ```

#### Protection Rules:
- **Required reviewers**: Senior SRE or incident response team only
- **Emergency access**: Can deploy from any branch during incidents
- **Audit logging**: All rollback deployments logged and tracked

#### Required Secrets:
- Same as production environment
- `ROLLBACK_NOTIFICATION_WEBHOOK`: Emergency notification channel

## Environment Security Best Practices

### Secrets Management
1. **Principle of Least Privilege**: Only add secrets that are absolutely necessary
2. **Regular Rotation**: Rotate API keys and credentials quarterly
3. **Environment Isolation**: Never share secrets between staging and production
4. **Audit Trail**: Monitor secret access and usage

### Access Control
1. **Team-Based Access**: Use GitHub teams for reviewer assignments
2. **Two-Person Rule**: Require multiple reviewers for production
3. **Branch Protection**: Restrict deployment branches for production
4. **Emergency Procedures**: Document emergency access procedures

### Monitoring and Alerting
1. **Deployment Notifications**: Configure Slack/Teams notifications
2. **Failed Deployment Alerts**: Monitor failed deployments
3. **Unauthorized Access**: Alert on unauthorized environment access
4. **Compliance Reporting**: Generate deployment audit reports

## Workflow Integration

### Environment Usage in Workflows

```yaml
# Staging deployment
deploy-staging:
  runs-on: ubuntu-latest
  environment:
    name: staging
    url: https://vcf-agent-staging.example.com
  steps:
    - name: Deploy to staging
      run: echo "Deploying to staging..."

# Production deployment with approval
deploy-production:
  runs-on: ubuntu-latest
  environment:
    name: production
    url: https://vcf-agent.example.com
  steps:
    - name: Deploy to production
      run: echo "Deploying to production..."

# Emergency rollback
rollback:
  runs-on: ubuntu-latest
  environment:
    name: production-rollback
    url: https://vcf-agent.example.com
  if: failure()
  steps:
    - name: Emergency rollback
      run: echo "Rolling back..."
```

### Environment Variables Access

```yaml
steps:
  - name: Deploy with environment variables
    env:
      ENVIRONMENT: ${{ vars.ENVIRONMENT }}
      LOG_LEVEL: ${{ vars.LOG_LEVEL }}
      DATABASE_URL: ${{ secrets.DATABASE_URL }}
    run: |
      echo "Deploying to $ENVIRONMENT"
      echo "Log level: $LOG_LEVEL"
      # Database URL is automatically available as secret
```

## Troubleshooting

### Common Issues

1. **Environment not found**
   - Verify environment name matches workflow exactly
   - Check repository permissions and access

2. **Approval required but no reviewers**
   - Add reviewers to environment configuration
   - Verify reviewer team membership

3. **Deployment branch restrictions**
   - Check branch protection rules
   - Verify deployment is from allowed branch/tag

4. **Secret access errors**
   - Verify secret exists in environment
   - Check secret name spelling and case sensitivity

### Linter Warnings

If your IDE or linter shows environment validation warnings:

1. **GitHub Actions linter**: May not recognize environments until they exist
2. **VS Code**: Install GitHub Actions extension for better validation
3. **Pre-commit hooks**: Configure to skip environment validation

### Emergency Procedures

If immediate deployment is required:

1. **Bypass environment protection** (Admin only):
   - Repository Settings > Environments > [environment] > Bypass protection rules
   - Document bypass reason in deployment notes

2. **Emergency rollback**:
   - Use `production-rollback` environment
   - Follow incident response procedures
   - Document rollback in incident log

## Validation Checklist

Before enabling production deployments:

- [ ] All three environments created (staging, production, production-rollback)
- [ ] Environment URLs configured correctly
- [ ] Required reviewers assigned and notified
- [ ] All required secrets added and tested
- [ ] Branch protection rules configured
- [ ] Notification webhooks tested
- [ ] Emergency procedures documented
- [ ] Team permissions verified
- [ ] Deployment workflow tested in staging

## References

- [GitHub Environments Documentation](https://docs.github.com/en/actions/deployment/targeting-different-environments/using-environments-for-deployment)
- [Environment Protection Rules](https://docs.github.com/en/actions/deployment/targeting-different-environments/using-environments-for-deployment#environment-protection-rules)
- [GitHub Actions Security](https://docs.github.com/en/actions/security-guides/security-hardening-for-github-actions)
- [Deployment Best Practices](https://docs.github.com/en/actions/deployment/about-deployments/deploying-with-github-actions)

---

**Next Steps**: After environment setup, test the deployment workflow in staging before enabling production deployments. 