{ pkgs, lib, config, inputs, ... }:

{
  # https://devenv.sh/basics/
  env.GREET = "Good day, happy coding";
  
  # Set R environment variables for better C++ integration
  env.R_LIBS_USER = "${config.env.DEVENV_STATE}/R";
  env.PKG_CONFIG_PATH = "${pkgs.pkg-config}/lib/pkgconfig";

  # https://devenv.sh/packages/
  packages = [ 
    pkgs.git

    # R development with C++ support
    pkgs.R
    pkgs.rPackages.Rcpp
    pkgs.rPackages.RcppEigen
    
    # VS Code R integration packages
    pkgs.rPackages.httpgd          # Plot viewer for VS Code
    pkgs.rPackages.languageserver  # R Language Server Protocol
    pkgs.rPackages.jsonlite        # JSON support for LSP
    pkgs.rPackages.renv            # Package management
 
    pkgs.rPackages.ggplot2         # Data visualization
    pkgs.rPackages.dplyr           # Data manipulation
    pkgs.rPackages.tidyr           # Data tidying
    
    # C++ development tools
    pkgs.gcc
    pkgs.gfortran
    pkgs.pkg-config
    pkgs.cmake
    pkgs.clang-tools_16
    
    
    # Linear algebra libraries (for RcppEigen)
    pkgs.eigen
    pkgs.blas
    pkgs.lapack
    
    # Additional R packages commonly used with Rcpp
    pkgs.rPackages.devtools
    pkgs.rPackages.testthat
    pkgs.rPackages.roxygen2
    pkgs.rPackages.knitr
  ];

  # https://devenv.sh/scripts/
  scripts.hello.exec = ''
    echo hello from $GREET
  '';
  
  scripts.r-setup.exec = ''
    echo "Setting up R environment..."
    mkdir -p $R_LIBS_USER
    echo "R library path: $R_LIBS_USER"
    echo "You can now use R with Rcpp and RcppEigen support!"
  '';

  scripts.test-rcpp.exec = ''
    echo "Testing Rcpp installation..."
    R --slave -e "library(Rcpp); cat('Rcpp version:', as.character(packageVersion('Rcpp')), '\n')"
    R --slave -e "library(RcppEigen); cat('RcppEigen version:', as.character(packageVersion('RcppEigen')), '\n')"
  '';

  scripts.setup-clangd.exec = ''
    echo "Setting up clangd configuration..."
    
    # Get actual R include paths dynamically
    R_INCLUDE=$(R --slave -e "cat(R.home('include'))")
    RCPP_INCLUDE=$(R --slave -e "cat(system.file('include', package='Rcpp'))")
    RCPPEIGEN_INCLUDE=$(R --slave -e "cat(system.file('include', package='RcppEigen'))")
    
    # Generate .clangd file with correct paths
    cat > .clangd <<EOF
CompileFlags:
  Add: [
    "-I./include",
    "-I${pkgs.eigen}/include/eigen3",
    "-I$R_INCLUDE",
    "-I$RCPP_INCLUDE", 
    "-I$RCPPEIGEN_INCLUDE"
  ]
EOF

    echo "✅ Generated .clangd with dynamic include paths:"
    cat .clangd
  '';

enterShell = ''
  hello
  git --version
  r-setup
  echo "C++ compiler: $(gcc --version | head -n1)"
  echo "R version: $(R --version | head -n1)"
  test-rcpp
  setup-clangd
  
  echo ""
  echo "🚀 R + C++ development environment is ready for VS Code!"
  echo "   - httpgd: Plot viewing"
  echo "   - languageserver: IntelliSense and code completion"
  echo "   - renv: Package management"
  echo "   - clangd: C++ LSP with correct R/Rcpp paths"
  echo ""
  echo "💡 If you still have issues, try running 'setup-clangd' manually"
  echo "   or restart your VS Code language server"
  echo ""
'';

  # https://devenv.sh/languages/
  languages.r.enable = true;
}