asm_Mol2SVG <- function(mol,
                       width = 500,
                       height = 500,
                       margin = 5,
                       zoom=2) {
  
  # Generate 2D coordinates
  generate.2d.coordinates(mol)
  
  # Create CDK depiction generator
  dg_class <- J("org.openscience.cdk.depict.DepictionGenerator")
  dg <- new(dg_class)
  
  # Set size
  dg <- dg$withSize(width, height)$withMargin(margin)$withZoom(zoom)
  
  # Generate depiction
  depiction <- dg$depict(mol)
  
  # Return SVG XML
  depiction$toSvgStr()
}



