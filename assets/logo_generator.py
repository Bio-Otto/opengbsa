#!/usr/bin/env python3
"""
OpenGBSA Logo Generator
Creates a simple and clean logo for the OpenGBSA package
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Circle, Rectangle
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import os

def create_svg_logo():
    """Create simple SVG logo for OpenGBSA"""
    
    # Create figure with transparent background
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 4)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # No background - transparent
    
    # Simple molecular representation - just a clean line
    x_line = np.linspace(1, 7, 30)
    y_line = 2 + 0.3 * np.sin(x_line * 2 * np.pi / 4)
    
    # Draw simple backbone
    ax.plot(x_line, y_line, color='#4CAF50', linewidth=3, alpha=0.8)
    
    # Add simple dots for residues
    for i in range(0, len(x_line), 6):
        if i < len(x_line):
            circle = Circle((x_line[i], y_line[i]), 0.08, 
                          facecolor='#2196F3', edgecolor='none')
            ax.add_patch(circle)
    
    # Simple text with black color
    ax.text(4, 3.2, 'OpenGBSA', fontsize=20, color='#000000', 
            weight='bold', ha='center', family='sans-serif')
    
    # Simple subtitle with black color
    ax.text(4, 0.8, 'MM/GBSA Analysis', fontsize=10, color='#333333', 
            ha='center', family='sans-serif')
    
    # Save as SVG
    plt.savefig('logo.svg', format='svg', dpi=300, bbox_inches='tight', 
                facecolor='none', edgecolor='none')
    plt.close()
    
    print("Simple SVG logo created: logo.svg")

def create_png_logo():
    """Create simple PNG logo for OpenGBSA"""
    
    # Create figure with transparent background
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # No background - transparent
    
    # Simple molecular representation
    x_line = np.linspace(1.5, 8.5, 40)
    y_line = 2.5 + 0.4 * np.sin(x_line * 2 * np.pi / 5)
    
    # Draw simple backbone
    ax.plot(x_line, y_line, color='#4CAF50', linewidth=4, alpha=0.9)
    
    # Add simple dots for residues
    for i in range(0, len(x_line), 8):
        if i < len(x_line):
            circle = Circle((x_line[i], y_line[i]), 0.12, 
                          facecolor='#2196F3', edgecolor='none')
            ax.add_patch(circle)
    
    # Simple text with black color
    ax.text(5, 3.8, 'OpenGBSA', fontsize=28, color='#000000', 
            weight='bold', ha='center', family='sans-serif')
    
    # Simple subtitle with black color
    ax.text(5, 1.2, 'MM/GBSA Analysis', fontsize=12, color='#333333', 
            ha='center', family='sans-serif')
    
    # Save as PNG
    plt.savefig('logo.png', format='png', dpi=300, bbox_inches='tight', 
                facecolor='none', edgecolor='none', transparent=True)
    plt.close()
    
    print("Simple PNG logo created: logo.png")

def create_simple_logo():
    """Create a very simple text-based logo"""
    
    # Create image with transparent background
    width, height = 800, 200
    img = Image.new('RGBA', (width, height), (0, 0, 0, 0))  # Transparent background
    draw = ImageDraw.Draw(img)
    
    try:
        font_large = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 48)
        font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 20)
    except:
        font_large = ImageFont.load_default()
        font_small = ImageFont.load_default()
    
    # Draw simple OpenGBSA text with black color
    text = "OpenGBSA"
    bbox = draw.textbbox((0, 0), text, font=font_large)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    
    x = (width - text_width) // 2
    y = (height - text_height) // 2 - 20
    
    draw.text((x, y), text, font=font_large, fill=(0, 0, 0, 255))  # Black text
    
    # Simple subtitle with dark gray color
    subtitle = "MM/GBSA Analysis"
    bbox_sub = draw.textbbox((0, 0), subtitle, font=font_small)
    sub_width = bbox_sub[2] - bbox_sub[0]
    
    x_sub = (width - sub_width) // 2
    y_sub = y + text_height + 10
    
    draw.text((x_sub, y_sub), subtitle, font=font_small, fill=(51, 51, 51, 255))  # Dark gray
    
    # Save
    img.save('logo_simple.png', 'PNG')
    print("Simple text logo created: logo_simple.png")

def create_github_banner():
    """Create a simple GitHub banner"""
    
    # Create image with GitHub banner dimensions
    width, height = 1280, 640
    img = Image.new('RGBA', (width, height), (13, 17, 23, 255))  # GitHub dark
    draw = ImageDraw.Draw(img)
    
    try:
        font_large = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 72)
        font_medium = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 32)
        font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 20)
    except:
        font_large = ImageFont.load_default()
        font_medium = ImageFont.load_default()
        font_small = ImageFont.load_default()
    
    # Draw simple OpenGBSA text with white color (for dark background)
    text = "OpenGBSA"
    bbox = draw.textbbox((0, 0), text, font=font_large)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    
    x = (width - text_width) // 2
    y = (height - text_height) // 2 - 60
    
    draw.text((x, y), text, font=font_large, fill=(255, 255, 255, 255))  # White text for dark background
    
    # Simple subtitle
    subtitle = "Molecular Mechanics / Generalized Born Surface Area Analysis"
    bbox_sub = draw.textbbox((0, 0), subtitle, font=font_medium)
    sub_width = bbox_sub[2] - bbox_sub[0]
    
    x_sub = (width - sub_width) // 2
    y_sub = y + text_height + 20
    
    draw.text((x_sub, y_sub), subtitle, font=font_medium, fill=(139, 148, 158, 255))
    
    # Simple feature list
    features = [
        "• MM/GBSA Energy Calculation",
        "• Per-Residue Decomposition", 
        "• Protein-Ligand Analysis",
        "• Advanced Reporting"
    ]
    
    y_features = y_sub + 80
    for i, feature in enumerate(features):
        draw.text((width//2 - 200, y_features + i * 40), feature, 
                 font=font_small, fill=(139, 148, 158, 255))
    
    # Save
    img.save('github_banner.png', 'PNG')
    print("Simple GitHub banner created: github_banner.png")

def main():
    """Generate all simple logo variants"""
    print("Generating simple OpenGBSA logos...")
    
    # Create different logo variants
    create_svg_logo()
    create_png_logo()
    create_simple_logo()
    create_github_banner()
    
    print("\nSimple logo generation complete!")
    print("Generated files:")
    print("  • logo.svg - Simple vector logo")
    print("  • logo.png - Simple high-resolution logo")
    print("  • logo_simple.png - Simple text logo")
    print("  • github_banner.png - Simple GitHub repository banner")
    
    print("\nUsage:")
    print("  • Add logo.png to your README.md")
    print("  • Use github_banner.png as your repository banner")
    print("  • Use logo.svg for web applications")

if __name__ == "__main__":
    main() 