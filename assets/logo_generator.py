#!/usr/bin/env python3
"""
OpenGBSA Logo Generator
Creates a professional and scientific logo for the OpenGBSA package
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import os

def create_svg_logo():
    """Create professional SVG logo for OpenGBSA"""
    
    # Create figure with professional dimensions
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Professional background with subtle gradient
    gradient = np.linspace(0, 1, 50)
    for i, alpha in enumerate(gradient):
        y_pos = i * 6 / 50
        color = (0.05, 0.1, 0.2, alpha * 0.15)
        rect = Rectangle((0, y_pos), 12, 6/50, facecolor=color, edgecolor='none')
        ax.add_patch(rect)
    
    # Professional molecular structure representation
    # Protein backbone with alpha helix representation
    x_helix = np.linspace(1.5, 10.5, 60)
    y_helix = 3 + 0.6 * np.sin(x_helix * 2 * np.pi / 3)
    
    # Draw professional backbone
    ax.plot(x_helix, y_helix, color='#2E7D32', linewidth=2.5, alpha=0.9)
    
    # Add professional residue representations
    for i in range(0, len(x_helix), 5):
        if i < len(x_helix):
            # Main residue
            circle = Circle((x_helix[i], y_helix[i]), 0.15, 
                          facecolor='#1976D2', edgecolor='#0D47A1', linewidth=1.5)
            ax.add_patch(circle)
            
            # Side chain representation
            if i % 10 == 0:
                x_side = x_helix[i] + 0.3 * np.cos(i * 0.5)
                y_side = y_helix[i] + 0.3 * np.sin(i * 0.5)
                side_circle = Circle((x_side, y_side), 0.08, 
                                   facecolor='#FF5722', edgecolor='#D84315', linewidth=1)
                ax.add_patch(side_circle)
    
    # Professional energy calculation symbols
    # Delta G symbol with scientific notation
    ax.text(3, 4.8, r'$\Delta G_{bind}$', fontsize=18, color='#FFC107', 
            weight='bold', ha='center', family='serif')
    
    # MM/GBSA text with professional styling
    ax.text(9, 4.8, 'MM/GBSA', fontsize=14, color='#FFFFFF', 
            weight='bold', ha='center', family='sans-serif')
    
    # Professional OpenGBSA text
    ax.text(6, 1.8, 'OpenGBSA', fontsize=26, color='#FFFFFF', 
            weight='bold', ha='center', family='sans-serif')
    
    # Professional subtitle
    ax.text(6, 1.2, 'Molecular Mechanics / Generalized Born Surface Area', 
            fontsize=9, color='#BDBDBD', ha='center', family='sans-serif')
    
    # Professional energy calculation visualization
    # Energy landscape representation
    x_energy = np.linspace(2, 10, 40)
    y_energy = 2.2 + 0.3 * np.sin(x_energy * 2 * np.pi / 4) + 0.1 * np.sin(x_energy * 4 * np.pi / 4)
    
    ax.plot(x_energy, y_energy, color='#FF9800', linewidth=2, alpha=0.8)
    
    # Add energy minima points
    for i in range(0, len(x_energy), 8):
        if i < len(x_energy):
            energy_point = Circle((x_energy[i], y_energy[i]), 0.06, 
                                facecolor='#FF9800', edgecolor='#F57C00', linewidth=1)
            ax.add_patch(energy_point)
    
    # Professional molecular interaction indicators
    for i in range(2):
        x_center = 4 + i * 4
        y_center = 1.5
        radius = 0.25
        
        # Interaction circle
        circle = Circle((x_center, y_center), radius, 
                       facecolor='none', edgecolor='#9C27B0', linewidth=1.5, alpha=0.7)
        ax.add_patch(circle)
        
        # Interaction center
        center = Circle((x_center, y_center), 0.05, 
                       facecolor='#9C27B0', edgecolor='none')
        ax.add_patch(center)
    
    # Save as SVG
    plt.savefig('logo.svg', format='svg', dpi=300, bbox_inches='tight', 
                facecolor='none', edgecolor='none')
    plt.close()
    
    print("Professional SVG logo created: logo.svg")

def create_png_logo():
    """Create professional PNG logo for OpenGBSA"""
    
    # Create figure with professional design
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Professional background with gradient
    gradient = np.linspace(0, 1, 60)
    for i, alpha in enumerate(gradient):
        y_pos = i * 8 / 60
        color = (0.03, 0.08, 0.15, alpha * 0.2)
        rect = Rectangle((0, y_pos), 14, 8/60, facecolor=color, edgecolor='none')
        ax.add_patch(rect)
    
    # Professional molecular structure
    # Alpha helix representation
    x_helix = np.linspace(2, 12, 70)
    y_helix = 4 + 0.8 * np.sin(x_helix * 2 * np.pi / 4)
    
    # Draw professional backbone
    ax.plot(x_helix, y_helix, color='#2E7D32', linewidth=3, alpha=0.9)
    
    # Add professional residue representations
    for i in range(0, len(x_helix), 6):
        if i < len(x_helix):
            # Main residue
            circle = Circle((x_helix[i], y_helix[i]), 0.2, 
                          facecolor='#1976D2', edgecolor='#0D47A1', linewidth=2)
            ax.add_patch(circle)
            
            # Side chain representation
            if i % 12 == 0:
                x_side = x_helix[i] + 0.4 * np.cos(i * 0.4)
                y_side = y_helix[i] + 0.4 * np.sin(i * 0.4)
                side_circle = Circle((x_side, y_side), 0.12, 
                                   facecolor='#FF5722', edgecolor='#D84315', linewidth=1.5)
                ax.add_patch(side_circle)
    
    # Professional energy calculation symbols
    # Delta G symbol
    ax.text(4, 6.5, r'$\Delta G_{bind}$', fontsize=24, color='#FFC107', 
            weight='bold', ha='center', family='serif')
    
    # MM/GBSA text
    ax.text(10, 6.5, 'MM/GBSA', fontsize=18, color='#FFFFFF', 
            weight='bold', ha='center', family='sans-serif')
    
    # Professional OpenGBSA text
    ax.text(7, 2.5, 'OpenGBSA', fontsize=32, color='#FFFFFF', 
            weight='bold', ha='center', family='sans-serif')
    
    # Professional subtitle
    ax.text(7, 1.8, 'Molecular Mechanics / Generalized Born Surface Area', 
            fontsize=11, color='#BDBDBD', ha='center', family='sans-serif')
    
    # Professional energy landscape
    x_energy = np.linspace(3, 11, 50)
    y_energy = 3.2 + 0.4 * np.sin(x_energy * 2 * np.pi / 5) + 0.15 * np.sin(x_energy * 3 * np.pi / 5)
    
    ax.plot(x_energy, y_energy, color='#FF9800', linewidth=2.5, alpha=0.8)
    
    # Add energy minima points
    for i in range(0, len(x_energy), 10):
        if i < len(x_energy):
            energy_point = Circle((x_energy[i], y_energy[i]), 0.08, 
                                facecolor='#FF9800', edgecolor='#F57C00', linewidth=1.5)
            ax.add_patch(energy_point)
    
    # Professional molecular interaction indicators
    for i in range(3):
        x_center = 4.5 + i * 2.5
        y_center = 2.2
        radius = 0.3
        
        # Interaction circle
        circle = Circle((x_center, y_center), radius, 
                       facecolor='none', edgecolor='#9C27B0', linewidth=2, alpha=0.7)
        ax.add_patch(circle)
        
        # Interaction center
        center = Circle((x_center, y_center), 0.06, 
                       facecolor='#9C27B0', edgecolor='none')
        ax.add_patch(center)
    
    # Save as PNG
    plt.savefig('logo.png', format='png', dpi=300, bbox_inches='tight', 
                facecolor='none', edgecolor='none', transparent=True)
    plt.close()
    
    print("Professional PNG logo created: logo.png")

def create_simple_logo():
    """Create a professional simple text-based logo"""
    
    # Create image with professional background
    width, height = 900, 250
    img = Image.new('RGBA', (width, height), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    
    try:
        font_large = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 56)
        font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 22)
    except:
        font_large = ImageFont.load_default()
        font_small = ImageFont.load_default()
    
    # Draw professional OpenGBSA text
    text = "OpenGBSA"
    bbox = draw.textbbox((0, 0), text, font=font_large)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    
    x = (width - text_width) // 2
    y = (height - text_height) // 2 - 25
    
    # Professional shadow effect
    draw.text((x+2, y+2), text, font=font_large, fill=(0, 0, 0, 80))
    draw.text((x, y), text, font=font_large, fill=(255, 255, 255, 255))
    
    # Professional subtitle
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
    """Create a professional GitHub banner"""
    
    # Create image with GitHub banner dimensions
    width, height = 1280, 640
    img = Image.new('RGBA', (width, height), (13, 17, 23, 255))
    draw = ImageDraw.Draw(img)
    
    try:
        font_large = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 80)
        font_medium = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 38)
        font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 26)
    except:
        font_large = ImageFont.load_default()
        font_medium = ImageFont.load_default()
        font_small = ImageFont.load_default()
    
    # Draw professional OpenGBSA text
    text = "OpenGBSA"
    bbox = draw.textbbox((0, 0), text, font=font_large)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    
    x = (width - text_width) // 2
    y = (height - text_height) // 2 - 60
    
    # Professional shadow effect
    draw.text((x+3, y+3), text, font=font_large, fill=(0, 0, 0, 100))
    draw.text((x, y), text, font=font_large, fill=(255, 255, 255, 255))
    
    # Professional subtitle
    subtitle = "Molecular Mechanics / Generalized Born Surface Area Analysis"
    bbox_sub = draw.textbbox((0, 0), subtitle, font=font_medium)
    sub_width = bbox_sub[2] - bbox_sub[0]
    
    x_sub = (width - sub_width) // 2
    y_sub = y + text_height + 25
    
    draw.text((x_sub, y_sub), subtitle, font=font_medium, fill=(139, 148, 158, 255))
    
    # Professional feature list
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