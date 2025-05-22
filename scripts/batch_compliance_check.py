#!/usr/bin/env python3
"""
Batch Compliance Checker for VCF/BCF Files
- Scans a directory for VCF/BCF files
- Optionally generates synthetic edge-case files
- Runs compliance checks (bcftools, GATK, or as configured)
- Outputs a markdown summary report
"""
import os
import sys
import argparse
from vcf_agent.validation import validate_vcf_file
import datetime

EDGE_CASES = [
    ("edgecase_missing_header.vcf", "#CHROM only, no meta headers", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00098\n22\t16051347\trs62224610\tG\tC\t0\tPASS\tAC=1;AF=0.5;AN=2\tGT:AD:DP:GD:GL:GQ:OG\t0|1:2,1:2:.:-3.95,-0.60,-3.69:32.52:./.\n"),
    ("edgecase_bad_info.vcf", "Malformed INFO field", "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00098\n22\t16051347\trs62224610\tG\tC\t0\tPASS\tAC=;AF=not_a_float;AN=2\tGT:AD:DP:GD:GL:GQ:OG\t0|1:2,1:2:.:-3.95,-0.60,-3.69:32.52:./.\n"),
    ("edgecase_nonstandard_chrom.vcf", "Non-standard chromosome name", "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00098\nchrUn\t16051347\trs62224610\tG\tC\t0\tPASS\tAC=1;AF=0.5;AN=2\tGT:AD:DP:GD:GL:GQ:OG\t0|1:2,1:2:.:-3.95,-0.60,-3.69:32.52:./.\n"),
    ("edgecase_symbolic_allele.vcf", "Symbolic ALT allele", "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00098\n22\t16051347\trs62224610\tG\t<DEL>\t0\tPASS\tAC=1;AF=0.5;AN=2\tGT:AD:DP:GD:GL:GQ:OG\t0|1:2,1:2:.:-3.95,-0.60,-3.69:32.52:./.\n"),
    ("edgecase_multiallelic.vcf", "Multi-allelic site", "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00098\n22\t16051347\trs62224610\tA\tG,T\t0\tPASS\tAC=1,2;AF=0.5,0.5;AN=2\tGT:AD:DP:GD:GL:GQ:OG\t1|2:2,1,1:2:.:-3.95,-0.60,-3.69:32.52:./.\n"),
    ("edgecase_missing_format.vcf", "Missing FORMAT column", "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tHG00098\n22\t16051347\trs62224610\tG\tC\t0\tPASS\tAC=1;AF=0.5;AN=2\t0|1:2,1:2:.:-3.95,-0.60,-3.69:32.52:./.\n"),
]

def generate_edgecases(directory):
    for fname, desc, content in EDGE_CASES:
        path = os.path.join(directory, fname)
        with open(path, "w") as f:
            f.write(content)
    print(f"Generated {len(EDGE_CASES)} edge-case VCF files in {directory}")

def is_vcf_or_bcf(filename):
    ext = filename.lower()
    return ext.endswith(".vcf") or ext.endswith(".vcf.gz") or ext.endswith(".bcf")

def generate_markdown_report(results, directory):
    lines = ["# Batch Compliance Check Report\n",
             f"Checked directory: `{directory}`\n",
             "| File | Tool | Valid | Error |",
             "|------|------|-------|-------|"]
    for fname, tool, valid, error in results:
        lines.append(f"| `{fname}` | `{tool}` | {'✅' if valid else '❌'} | {error.replace('|', ' ')} |")
    return "\n".join(lines)

def generate_color_markdown_report(results, directory):
    # Uses HTML <span> for color, which works in GitHub/HTML renderers
    lines = ["# Batch Compliance Check Report\n",
             f"Checked directory: `{directory}`\n",
             "| File | Tool | Valid | Error |",
             "|------|------|-------|-------|"]
    for fname, tool, valid, error in results:
        if valid:
            valid_str = '<span style="color:green;font-weight:bold">✅ PASS</span>'
        else:
            valid_str = '<span style="color:red;font-weight:bold">❌ FAIL</span>'
        lines.append(f"| `{fname}` | `{tool}` | {valid_str} | {error.replace('|', ' ')} |")
    # Add a legend
    lines.append("\n<sub>✅ <span style='color:green'>PASS</span>, ❌ <span style='color:red'>FAIL</span></sub>")
    return "\n".join(lines)

def generate_html_report(results, directory):
    now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    html = [
        "<!DOCTYPE html>",
        "<html lang='en'>",
        "<head>",
        "<meta charset='UTF-8'>",
        "<title>Batch Compliance Check Report</title>",
        "<style>",
        "body { font-family: Arial, sans-serif; background: #f9f9f9; color: #222; }",
        "h1 { color: #2c3e50; }",
        "table { border-collapse: collapse; width: 100%; background: #fff; }",
        "th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }",
        "th { background: #f2f2f2; }",
        ".pass { color: #27ae60; font-weight: bold; }",
        ".fail { color: #c0392b; font-weight: bold; }",
        ".legend { margin-top: 1em; font-size: 0.95em; }",
        "</style>",
        "</head>",
        "<body>",
        f"<h1>Batch Compliance Check Report</h1>",
        f"<p><b>Checked directory:</b> <code>{directory}</code><br><b>Generated:</b> {now}</p>",
        "<table>",
        "<tr><th>File</th><th>Tool</th><th>Valid</th><th>Error</th></tr>"
    ]
    for fname, tool, valid, error in results:
        if valid:
            valid_str = '<span class="pass">✅ PASS</span>'
        else:
            valid_str = '<span class="fail">❌ FAIL</span>'
        html.append(f"<tr><td><code>{fname}</code></td><td><code>{tool}</code></td><td>{valid_str}</td><td>{error.replace('|', ' ')}</td></tr>")
    html.extend([
        "</table>",
        "<div class='legend'>✅ <span class='pass'>PASS</span>, ❌ <span class='fail'>FAIL</span></div>",
        "</body>",
        "</html>"
    ])
    return "\n".join(html)

def main():
    parser = argparse.ArgumentParser(description="Batch VCF/BCF Compliance Checker")
    parser.add_argument("-d", "--directory", default="sample_data", help="Directory to scan for VCF/BCF files")
    parser.add_argument("-t", "--tool", default=None, help="Compliance tool to use (bcftools, gatk, or as configured)")
    parser.add_argument("-o", "--output", default=None, help="Output file (default: stdout, extension based on format)")
    parser.add_argument("--generate-edgecases", action="store_true", help="Generate synthetic edge-case VCF files in the directory")
    parser.add_argument("--notify", action="store_true", help="Save report to notifications/ folder (default filename: compliance_report.md or .html)")
    parser.add_argument("--format", choices=["markdown", "color-markdown", "html"], default="markdown", help="Output format: markdown, color-markdown, or html")
    args = parser.parse_args()

    if args.generate_edgecases:
        generate_edgecases(args.directory)

    files = [f for f in os.listdir(args.directory) if is_vcf_or_bcf(f)]
    if not files:
        print(f"No VCF/BCF files found in {args.directory}")
        return

    results = []
    for fname in files:
        path = os.path.join(args.directory, fname)
        try:
            is_valid, error = validate_vcf_file(path, tool=args.tool)
        except Exception as e:
            is_valid, error = False, f"Exception: {e}"
        results.append((fname, args.tool or "(default)", is_valid, error or ""))

    # Select report generator
    if args.format == "markdown":
        report = generate_markdown_report(results, args.directory)
        ext = ".md"
    elif args.format == "color-markdown":
        report = generate_color_markdown_report(results, args.directory)
        ext = ".md"
    elif args.format == "html":
        report = generate_html_report(results, args.directory)
        ext = ".html"
    else:
        raise ValueError("Unknown format")

    # Output logic
    if args.notify:
        notif_dir = "notifications"
        os.makedirs(notif_dir, exist_ok=True)
        base = args.output or f"compliance_report{ext}"
        notif_path = os.path.join(notif_dir, base)
        with open(notif_path, "w") as f:
            f.write(report)
        print(f"[NOTIFY] Compliance report saved to {notif_path}")
    elif args.output:
        with open(args.output, "w") as f:
            f.write(report)
        print(f"Report written to {args.output}")
    else:
        print(report)

if __name__ == "__main__":
    main() 