# How to contribute
I welcome contributions to this project! Whether it's bug fixes, new features, or documentation improvements, your help is appreciated. Here's how you can contribute:
1. **Fork the Repository**: Click the "Fork" button at the top right of this page to create your own copy of the repository.
2. **Clone Your Fork**: Clone your forked repository to your local machine using:
   ```bash
   git clone https://github.com/yourusername/repo-name.git
    cd repo-name
    ```
3. **Create a New Branch**: Create a new branch for your feature or bug fix:
    ```bash
    git checkout -b feature-or-bugfix-name
    ```
4. **Make Your Changes**: Implement your changes in the codebase. Ensure that your code follows the existing style and conventions.
5. **Test Your Changes**: Run existing tests and add new tests if necessary to ensure your changes work as expected.
6. **Commit Your Changes**: Commit your changes with a descriptive message:
    ```bash
    git add .
    git commit -m "Description of your changes"
    ```
7. **Push to Your Fork**: Push your changes to your forked repository:
    ```bash
    git push origin feature-or-bugfix-name
    ```
8. **Create a Pull Request**: Go to the original repository and create a pull request from your forked branch. Provide a clear description of the changes you made and any relevant context.
9. **Address Feedback**: Be responsive to any feedback or requests for changes from the maintainers.

# Tips&Tricks

- The philosophy of this project is to have modular components that can be easily combined to create custom MCMC schemes for Bayesian nonparametric models. My personal guidelines is to follow the **KISS** principle (Keep It Simple, Stupid) when adding new features or making changes. This means avoiding unnecessary complexity and focusing on clear, maintainable code.

- When adding new features, consider how they can be integrated into the existing modular structure. For example, if you're adding a new likelihood model, ensure it adheres to the `Likelihood` abstract class interface so it can be easily swapped with other likelihoods.

- Important is to keep the code well-documented, especially when introducing new components. Clear comments and documentation help others (and your future self) understand the purpose and functionality of the code.

- Once added a new feature, remember to extend the R bindings in `src/bindings.cpp` to expose the new functionality to R users.

