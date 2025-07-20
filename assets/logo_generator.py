#!/usr/bin/env python3
"""
OpenGBSA Logo Generator
Creates a professional logo for the OpenGBSA package
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import os

def create_svg_logo():
    """Create SVG logo for OpenGBSA"""
    
    # Create figure with transparent background
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 6)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Background gradient
    gradient = np.linspace(0, 1, 100)
    for i, alpha in enumerate(gradient):
        y_pos = i * 6 / 100
        color = (0.1, 0.2, 0.4, alpha * 0.3)  # Dark blue gradient
        rect = Rectangle((0, y_pos), 8, 6/100, facecolor=color, edgecolor='none')
        ax.add_patch(rect)
    
    # Main molecule representation (simplified DNA/protein structure)
    # Helix backbone
    x_helix = np.linspace(1, 7, 50)
    y_helix = 3 + 0.8 * np.sin(x_helix * 2 * np.pi / 3)
    
    # Draw helix backbone
    ax.plot(x_helix, y_helix, color='#4CAF50', linewidth=3, alpha=0.8)
    
    # Add base pairs (simplified as circles)
    for i in range(0, len(x_helix), 3):
        if i < len(x_helix):
            # Base pair 1
            circle1 = Circle((x_helix[i], y_helix[i] + 0.3), 0.15, 
                           facecolor='#2196F3', edgecolor='#1976D2', linewidth=2)
            ax.add_patch(circle1)
            
            # Base pair 2
            circle2 = Circle((x_helix[i], y_helix[i] - 0.3), 0.15, 
                           facecolor='#FF5722', edgecolor='#D84315', linewidth=2)
            ax.add_patch(circle2)
    
    # Add energy calculation symbols
    # Delta G symbol
    ax.text(2, 4.5, r'$\Delta G$', fontsize=24, color='#FFC107', 
            weight='bold', ha='center')
    
    # MM/GBSA text
    ax.text(6, 4.5, 'MM/GBSA', fontsize=16, color='#FFFFFF', 
            weight='bold', ha='center')
    
    # OpenGBSA text
    ax.text(4, 1.2, 'OpenGBSA', fontsize=28, color='#FFFFFF', 
            weight='bold', ha='center', family='sans-serif')
    
    # Subtitle
    ax.text(4, 0.6, 'Molecular Mechanics / Generalized Born Surface Area', 
            fontsize=10, color='#BDBDBD', ha='center', family='sans-serif')
    
    # Add some molecular energy lines
    for i in range(3):
        x_start = 1.5 + i * 2
        y_start = 2.5
        x_end = x_start + 1
        y_end = y_start + 0.3 * np.sin(i * np.pi / 2)
        
        ax.plot([x_start, x_end], [y_start, y_end], 
                color='#FF9800', linewidth=2, alpha=0.7)
    
    # Save as SVG
    plt.savefig('logo.svg', format='svg', dpi=300, bbox_inches='tight', 
                facecolor='none', edgecolor='none')
    plt.close()
    
    print("SVG logo created: logo.svg")

def create_png_logo():
    """Create PNG logo for OpenGBSA"""
    
    # Create figure with transparent background
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Background with gradient
    gradient = np.linspace(0, 1, 100)
    for i, alpha in enumerate(gradient):
        y_pos = i * 8 / 100
        color = (0.05, 0.1, 0.2, alpha * 0.4)  # Darker blue gradient
        rect = Rectangle((0, y_pos), 10, 8/100, facecolor=color, edgecolor='none')
        ax.add_patch(rect)
    
    # Main molecule representation
    # DNA/protein helix
    x_helix = np.linspace(1.5, 8.5, 60)
    y_helix = 4 + 1.0 * np.sin(x_helix * 2 * np.pi / 4)
    
    # Draw helix backbone
    ax.plot(x_helix, y_helix, color='#4CAF50', linewidth=4, alpha=0.9)
    
    # Add base pairs
    for i in range(0, len(x_helix), 4):
        if i < len(x_helix):
            # Base pair 1
            circle1 = Circle((x_helix[i], y_helix[i] + 0.4), 0.2, 
                           facecolor='#2196F3', edgecolor='#1976D2', linewidth=3)
            ax.add_patch(circle1)
            
            # Base pair 2
            circle2 = Circle((x_helix[i], y_helix[i] - 0.4), 0.2, 
                           facecolor='#FF5722', edgecolor='#D84315', linewidth=3)
            ax.add_patch(circle2)
    
    # Energy calculation symbols
    # Delta G symbol
    ax.text(3, 6.5, r'$\Delta G$', fontsize=32, color='#FFC107', 
            weight='bold', ha='center')
    
    # MM/GBSA text
    ax.text(7, 6.5, 'MM/GBSA', fontsize=20, color='#FFFFFF', 
            weight='bold', ha='center')
    
    # OpenGBSA text
    ax.text(5, 1.5, 'OpenGBSA', fontsize=36, color='#FFFFFF', 
            weight='bold', ha='center', family='sans-serif')
    
    # Subtitle
    ax.text(5, 0.8, 'Molecular Mechanics / Generalized Born Surface Area', 
            fontsize=12, color='#BDBDBD', ha='center', family='sans-serif')
    
    # Add energy calculation lines
    for i in range(4):
        x_start = 2 + i * 1.5
        y_start = 3.2
        x_end = x_start + 1.2
        y_end = y_start + 0.4 * np.sin(i * np.pi / 2)
        
        ax.plot([x_start, x_end], [y_start, y_end], 
                color='#FF9800', linewidth=3, alpha=0.8)
    
    # Add some molecular interaction lines
    for i in range(3):
        x_center = 3 + i * 2
        y_center = 2.5
        radius = 0.3
        
        circle = Circle((x_center, y_center), radius, 
                       facecolor='none', edgecolor='#9C27B0', linewidth=2, alpha=0.6)
        ax.add_patch(circle)
    
    # Save as PNG
    plt.savefig('logo.png', format='png', dpi=300, bbox_inches='tight', 
                facecolor='none', edgecolor='none', transparent=True)
    plt.close()
    
    print("PNG logo created: logo.png")

def create_simple_logo():
    """Create a simple text-based logo"""
    
    # Create image with transparent background
    width, height = 800, 200
    img = Image.new('RGBA', (width, height), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    
    try:
        # Try to use a nice font
        font_large = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 48)
        font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 20)
    except:
        # Fallback to default font
        font_large = ImageFont.load_default()
        font_small = ImageFont.load_default()
    
    # Draw OpenGBSA text
    text = "OpenGBSA"
    bbox = draw.textbbox((0, 0), text, font=font_large)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    
    x = (width - text_width) // 2
    y = (height - text_height) // 2 - 20
    
    # Draw text with shadow effect
    draw.text((x+2, y+2), text, font=font_large, fill=(0, 0, 0, 100))
    draw.text((x, y), text, font=font_large, fill=(255, 255, 255, 255))
    
    # Draw subtitle
    subtitle = "Molecular Mechanics / Generalized Born Surface Area"
    bbox_sub = draw.textbbox((0, 0), subtitle, font=font_small)
    sub_width = bbox_sub[2] - bbox_sub[0]
    
    x_sub = (width - sub_width) // 2
    y_sub = y + text_height + 10
    
    draw.text((x_sub, y_sub), subtitle, font=font_small, fill=(189, 189, 189, 255))
    
    # Save
    img.save('logo_simple.png', 'PNG')
    print("Simple logo created: logo_simple.png")

def create_github_banner():
    """Create a GitHub banner for the repository"""
    
    # Create image with GitHub banner dimensions
    width, height = 1280, 640
    img = Image.new('RGBA', (width, height), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    
    try:
        font_large = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 72)
        font_medium = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 36)
        font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 24)
    except:
        font_large = ImageFont.load_default()
        font_medium = ImageFont.load_default()
        font_small = ImageFont.load_default()
    
    # Draw gradient background
    for y in range(height):
        alpha = int(255 * (1 - y / height) * 0.3)
        color = (13, 17, 23, alpha)  # GitHub dark theme
        draw.line([(0, y), (width, y)], fill=color)
    
    # Draw OpenGBSA text
    text = "OpenGBSA"
    bbox = draw.textbbox((0, 0), text, font=font_large)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    
    x = (width - text_width) // 2
    y = (height - text_height) // 2 - 50
    
    # Draw text with glow effect
    for offset in range(5, 0, -1):
        alpha = 255 // (offset + 1)
        draw.text((x+offset, y+offset), text, font=font_large, fill=(255, 255, 255, alpha))
    
    draw.text((x, y), text, font=font_large, fill=(255, 255, 255, 255))
    
    # Draw subtitle
    subtitle = "Molecular Mechanics / Generalized Born Surface Area Analysis"
    bbox_sub = draw.textbbox((0, 0), subtitle, font=font_medium)
    sub_width = bbox_sub[2] - bbox_sub[0]
    
    x_sub = (width - sub_width) // 2
    y_sub = y + text_height + 20
    
    draw.text((x_sub, y_sub), subtitle, font=font_medium, fill=(139, 148, 158, 255))
    
    # Draw features
    features = [
        "• Complete MM/GBSA Energy Calculation",
        "• Per-Residue Energy Decomposition", 
        "• Protein-Ligand Interaction Analysis",
        "• Advanced Visualization & Reporting"
    ]
    
    y_features = y_sub + 80
    for i, feature in enumerate(features):
        draw.text((width//2 - 300, y_features + i * 40), feature, 
                 font=font_small, fill=(139, 148, 158, 255))
    
    # Save
    img.save('github_banner.png', 'PNG')
    print("GitHub banner created: github_banner.png")

def main():
    """Generate all logo variants"""
    print("Generating OpenGBSA logos...")
    
    # Create different logo variants
    create_svg_logo()
    create_png_logo()
    create_simple_logo()
    create_github_banner()
    
    print("\nLogo generation complete!")
    print("Generated files:")
    print("  • logo.svg - Vector logo for web")
    print("  • logo.png - High-resolution logo")
    print("  • logo_simple.png - Simple text logo")
    print("  • github_banner.png - GitHub repository banner")
    
    print("\nUsage:")
    print("  • Add logo.png to your README.md")
    print("  • Use github_banner.png as your repository banner")
    print("  • Use logo.svg for web applications")

if __name__ == "__main__":
    main() 