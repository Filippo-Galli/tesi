{ pkgs, lib, config, inputs, ... }:

{
  # https://devenv.sh/basics/
  env.GREET = "Good day, happy coding";
  
  # Set R environment variables for better C++ integration
  env.R_LIBS_USER = "${config.env.DEVENV_STATE}/R";
  env.PKG_CONFIG_PATH = "${pkgs.pkg-config}/lib/pkgconfig";
  
  # Ensure R can find Nix-installed packages
  env.R_LIBS_SITE = "${pkgs.R}/library";

  # https://devenv.sh/packages/
  packages = [ 
    pkgs.git

    # R development with C++ support
    pkgs.R
    pkgs.rPackages.Rcpp
    pkgs.rPackages.RcppEigen
    
    # VS Code R integration packages (remove duplicates)
    pkgs.rPackages.httpgd          # Plot viewer for VS Code
    pkgs.rPackages.languageserver  # R Language Server Protocol (keep only once)
    pkgs.rPackages.jsonlite        # JSON support for LSP
    pkgs.rPackages.renv            # Package management
 
    # Data analysis packages
    pkgs.rPackages.ggplot2         # Data visualization
    pkgs.rPackages.dplyr           # Data manipulation
    pkgs.rPackages.tidyr           # Data tidying

    # MCMC and clustering tools
    pkgs.rPackages.spam
    pkgs.rPackages.fields
    pkgs.rPackages.viridisLite
    pkgs.rPackages.RColorBrewer
    pkgs.rPackages.pheatmap
    pkgs.rPackages.mcclust
    pkgs.rPackages.salso           # SALSO clustering
    pkgs.rPackages.mclust          # For ARI calculation
    
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
    
    # Additional R development packages
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
    echo "R site library: $R_LIBS_SITE"
    echo "You can now use R with Rcpp, RcppEigen, and SALSO support!"
  '';

  scripts.test-rcpp.exec = ''
    echo "Testing R package installations..."
    R --slave -e "library(Rcpp); cat('Rcpp version:', as.character(packageVersion('Rcpp')), '\n')"
    R --slave -e "library(RcppEigen); cat('RcppEigen version:', as.character(packageVersion('RcppEigen')), '\n')"
    R --slave -e "library(salso); cat('SALSO version:', as.character(packageVersion('salso')), '\n')"
    R --slave -e "library(languageserver); cat('languageserver version:', as.character(packageVersion('languageserver')), '\n')"
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

  scripts.check-r-packages.exec = ''
    echo "Checking all required R packages..."
    R --slave -e "
      required_packages <- c('salso', 'pheatmap', 'mclust', 'mcclust', 'languageserver', 'httpgd')
      for (pkg in required_packages) {
        if (requireNamespace(pkg, quietly = TRUE)) {
          cat('✅', pkg, 'version:', as.character(packageVersion(pkg)), '\n')
        } else {
          cat('❌', pkg, 'not found\n')
        }
      }
    "
  '';

  scripts.clean_o_files.exec = ''
    echo "Cleaning up .o files..."
    find . -name "*.o" -type f -delete
    echo "✅ Removed all .o files"
  '';


enterShell = ''
  setup-clangd
  
  echo ""
  echo "🚀 R + C++ development environment is ready for VS Code!"
  echo "   - httpgd: Plot viewing"
  echo "   - languageserver: IntelliSense and code completion"
  echo "   - salso: Modern clustering analysis"
  echo "   - clangd: C++ LSP with correct R/Rcpp paths"
  echo ""
  echo "💡 If VS Code asks to install languageserver, click 'No' - it's already available via Nix"
  echo "   or restart your VS Code language server"
  echo "💡 If you see 'RcppEigen not found', run 'devenv r-setup' to set up the environment"
  echo "💡 to clean C++ objects file .o type clean_o_files"
  echo ""
'';

  # https://devenv.sh/languages/
  languages.r.enable = true;
}