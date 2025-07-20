#!/usr/bin/env python3
"""
OpenGBSA Logo Generator
Creates a professional logo for the OpenGBSA package
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle, Ellipse
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import os

def create_svg_logo():
    """Create professional SVG logo for OpenGBSA"""
    
    # Create figure with clean background
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Clean background with subtle gradient
    gradient = np.linspace(0, 1, 50)
    for i, alpha in enumerate(gradient):
        y_pos = i * 6 / 50
        color = (0.05, 0.1, 0.2, alpha * 0.2)  # Subtle dark blue
        rect = Rectangle((0, y_pos), 10, 6/50, facecolor=color, edgecolor='none')
        ax.add_patch(rect)
    
    # Professional molecular structure representation
    # Clean protein backbone
    x_backbone = np.linspace(1.5, 8.5, 40)
    y_backbone = 3 + 0.6 * np.sin(x_backbone * 2 * np.pi / 4)
    
    # Draw clean backbone
    ax.plot(x_backbone, y_backbone, color='#4CAF50', linewidth=2.5, alpha=0.8)
    
    # Add clean residue representations
    for i in range(0, len(x_backbone), 4):
        if i < len(x_backbone):
            # Residue 1
            circle1 = Circle((x_backbone[i], y_backbone[i] + 0.25), 0.12, 
                           facecolor='#2196F3', edgecolor='#1976D2', linewidth=1.5)
            ax.add_patch(circle1)
            
            # Residue 2
            circle2 = Circle((x_backbone[i], y_backbone[i] - 0.25), 0.12, 
                           facecolor='#FF5722', edgecolor='#D84315', linewidth=1.5)
            ax.add_patch(circle2)
    
    # Professional energy calculation symbols
    # Clean Delta G symbol
    ax.text(3, 4.8, r'$\Delta G$', fontsize=20, color='#FFC107', 
            weight='bold', ha='center', family='serif')
    
    # Clean MM/GBSA text
    ax.text(7, 4.8, 'MM/GBSA', fontsize=14, color='#FFFFFF', 
            weight='bold', ha='center', family='sans-serif')
    
    # Professional OpenGBSA text with clean layout
    ax.text(5, 1.8, 'OpenGBSA', fontsize=24, color='#FFFFFF', 
            weight='bold', ha='center', family='sans-serif')
    
    # Clean subtitle
    ax.text(5, 1.2, 'Molecular Mechanics / Generalized Born Surface Area', 
            fontsize=9, color='#BDBDBD', ha='center', family='sans-serif')
    
    # Professional energy calculation lines
    for i in range(3):
        x_start = 2.5 + i * 2.5
        y_start = 2.8
        x_end = x_start + 1.5
        y_end = y_start + 0.2 * np.sin(i * np.pi / 2)
        
        ax.plot([x_start, x_end], [y_start, y_end], 
                color='#FF9800', linewidth=2, alpha=0.7)
    
    # Add professional molecular interaction indicators
    for i in range(2):
        x_center = 3.5 + i * 3
        y_center = 2.2
        radius = 0.25
        
        circle = Circle((x_center, y_center), radius, 
                       facecolor='none', edgecolor='#9C27B0', linewidth=1.5, alpha=0.6)
        ax.add_patch(circle)
    
    # Save as SVG
    plt.savefig('logo.svg', format='svg', dpi=300, bbox_inches='tight', 
                facecolor='none', edgecolor='none')
    plt.close()
    
    print("Professional SVG logo created: logo.svg")

def create_png_logo():
    """Create professional PNG logo for OpenGBSA"""
    
    # Create figure with clean design
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Professional background with clean gradient
    gradient = np.linspace(0, 1, 60)
    for i, alpha in enumerate(gradient):
        y_pos = i * 8 / 60
        color = (0.03, 0.08, 0.15, alpha * 0.25)  # Clean dark blue
        rect = Rectangle((0, y_pos), 12, 8/60, facecolor=color, edgecolor='none')
        ax.add_patch(rect)
    
    # Professional molecular structure
    # Clean protein backbone
    x_backbone = np.linspace(2, 10, 50)
    y_backbone = 4 + 0.8 * np.sin(x_backbone * 2 * np.pi / 5)
    
    # Draw clean backbone
    ax.plot(x_backbone, y_backbone, color='#4CAF50', linewidth=3, alpha=0.9)
    
    # Add clean residue representations
    for i in range(0, len(x_backbone), 5):
        if i < len(x_backbone):
            # Residue 1
            circle1 = Circle((x_backbone[i], y_backbone[i] + 0.35), 0.18, 
                           facecolor='#2196F3', edgecolor='#1976D2', linewidth=2)
            ax.add_patch(circle1)
            
            # Residue 2
            circle2 = Circle((x_backbone[i], y_backbone[i] - 0.35), 0.18, 
                           facecolor='#FF5722', edgecolor='#D84315', linewidth=2)
            ax.add_patch(circle2)
    
    # Professional energy calculation symbols
    # Clean Delta G symbol
    ax.text(4, 6.5, r'$\Delta G$', fontsize=28, color='#FFC107', 
            weight='bold', ha='center', family='serif')
    
    # Clean MM/GBSA text
    ax.text(8, 6.5, 'MM/GBSA', fontsize=18, color='#FFFFFF', 
            weight='bold', ha='center', family='sans-serif')
    
    # Professional OpenGBSA text
    ax.text(6, 2.2, 'OpenGBSA', fontsize=32, color='#FFFFFF', 
            weight='bold', ha='center', family='sans-serif')
    
    # Clean subtitle
    ax.text(6, 1.5, 'Molecular Mechanics / Generalized Born Surface Area', 
            fontsize=11, color='#BDBDBD', ha='center', family='sans-serif')
    
    # Professional energy calculation lines
    for i in range(4):
        x_start = 3 + i * 2
        y_start = 3.5
        x_end = x_start + 1.8
        y_end = y_start + 0.3 * np.sin(i * np.pi / 2)
        
        ax.plot([x_start, x_end], [y_start, y_end], 
                color='#FF9800', linewidth=2.5, alpha=0.8)
    
    # Add professional molecular interaction indicators
    for i in range(3):
        x_center = 4 + i * 2.5
        y_center = 2.8
        radius = 0.3
        
        circle = Circle((x_center, y_center), radius, 
                       facecolor='none', edgecolor='#9C27B0', linewidth=2, alpha=0.6)
        ax.add_patch(circle)
    
    # Save as PNG
    plt.savefig('logo.png', format='png', dpi=300, bbox_inches='tight', 
                facecolor='none', edgecolor='none', transparent=True)
    plt.close()
    
    print("Professional PNG logo created: logo.png")

def create_simple_logo():
    """Create a clean and professional simple text-based logo"""
    
    # Create image with clean background
    width, height = 900, 250
    img = Image.new('RGBA', (width, height), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    
    try:
        # Use professional fonts
        font_large = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 56)
        font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 22)
    except:
        font_large = ImageFont.load_default()
        font_small = ImageFont.load_default()
    
    # Draw clean OpenGBSA text
    text = "OpenGBSA"
    bbox = draw.textbbox((0, 0), text, font=font_large)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    
    x = (width - text_width) // 2
    y = (height - text_height) // 2 - 25
    
    # Clean shadow effect
    draw.text((x+2, y+2), text, font=font_large, fill=(0, 0, 0, 80))
    draw.text((x, y), text, font=font_large, fill=(255, 255, 255, 255))
    
    # Clean subtitle
    subtitle = "Molecular Mechanics / Generalized Born Surface Area"
    bbox_sub = draw.textbbox((0, 0), subtitle, font=font_small)
    sub_width = bbox_sub[2] - bbox_sub[0]
    
    x_sub = (width - sub_width) // 2
    y_sub = y + text_height + 15
    
    draw.text((x_sub, y_sub), subtitle, font=font_small, fill=(189, 189, 189, 255))
    
    # Save
    img.save('logo_simple.png', 'PNG')
    print("Professional simple logo created: logo_simple.png")

def create_github_banner():
    """Create a clean and professional GitHub banner"""
    
    # Create image with GitHub banner dimensions
    width, height = 1280, 640
    img = Image.new('RGBA', (width, height), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    
    try:
        font_large = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 80)
        font_medium = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 38)
        font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 26)
    except:
        font_large = ImageFont.load_default()
        font_medium = ImageFont.load_default()
        font_small = ImageFont.load_default()
    
    # Clean gradient background
    for y in range(height):
        alpha = int(255 * (1 - y / height) * 0.25)
        color = (13, 17, 23, alpha)  # Clean GitHub dark theme
        draw.line([(0, y), (width, y)], fill=color)
    
    # Draw clean OpenGBSA text
    text = "OpenGBSA"
    bbox = draw.textbbox((0, 0), text, font=font_large)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    
    x = (width - text_width) // 2
    y = (height - text_height) // 2 - 60
    
    # Clean shadow effect
    draw.text((x+3, y+3), text, font=font_large, fill=(0, 0, 0, 100))
    draw.text((x, y), text, font=font_large, fill=(255, 255, 255, 255))
    
    # Clean subtitle
    subtitle = "Molecular Mechanics / Generalized Born Surface Area Analysis"
    bbox_sub = draw.textbbox((0, 0), subtitle, font=font_medium)
    sub_width = bbox_sub[2] - bbox_sub[0]
    
    x_sub = (width - sub_width) // 2
    y_sub = y + text_height + 25
    
    draw.text((x_sub, y_sub), subtitle, font=font_medium, fill=(139, 148, 158, 255))
    
    # Clean feature list
    features = [
        "• Complete MM/GBSA Energy Calculation",
        "• Per-Residue Energy Decomposition", 
        "• Protein-Ligand Interaction Analysis",
        "• Advanced Visualization & Reporting"
    ]
    
    y_features = y_sub + 90
    for i, feature in enumerate(features):
        draw.text((width//2 - 320, y_features + i * 45), feature, 
                 font=font_small, fill=(139, 148, 158, 255))
    
    # Save
    img.save('github_banner.png', 'PNG')
    print("Professional GitHub banner created: github_banner.png")

def main():
    """Generate all professional logo variants"""
    print("Generating professional OpenGBSA logos...")
    
    # Create different logo variants
    create_svg_logo()
    create_png_logo()
    create_simple_logo()
    create_github_banner()
    
    print("\nProfessional logo generation complete!")
    print("Generated files:")
    print("  • logo.svg - Professional vector logo")
    print("  • logo.png - Professional high-resolution logo")
    print("  • logo_simple.png - Professional simple text logo")
    print("  • github_banner.png - Professional GitHub repository banner")
    
    print("\nUsage:")
    print("  • Add logo.png to your README.md")
    print("  • Use github_banner.png as your repository banner")
    print("  • Use logo.svg for web applications")

if __name__ == "__main__":
    main() 