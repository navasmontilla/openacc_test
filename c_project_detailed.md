# C Project: Macros via CMake and Python Build Automation

This project demonstrates how to:
- Structure a C project with headers and source files
- Use `CMake` to compile the code
- Define macros (e.g., `EQUATION_SYSTEM`, `USE_ADVANCED`) with default values in a header file
- Override macros using CMake flags
- Use a Python script to configure, compile, and run the project

---

## 🗂 Project Structure
```
main/
├── CMakeLists.txt
├── lib/
│   ├── definitions.h
│   ├── math_utils.c
│   └── math_utils.h
├── main.c
└── run.py
```

---

## 📄 `main/CMakeLists.txt`

Defines the build process and includes source and header files.
```cmake
cmake_minimum_required(VERSION 3.10)
project(MyCProject C)

# Include the lib directory for headers
include_directories(${PROJECT_SOURCE_DIR}/lib)

# Add the executable
add_executable(my_program main.c lib/math_utils.c)
```

---

## 📄 `main/lib/definitions.h`

Provides default values for macros, ensuring compilation even if CMake flags aren't set.
```c
#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// Default values for macros (overridable by -D flags)
#ifndef EQUATION_SYSTEM
#define EQUATION_SYSTEM 1
#endif

#ifndef USE_ADVANCED
#define USE_ADVANCED 0
#endif

#endif // DEFINITIONS_H
```

---

## 📄 `main/lib/math_utils.h`

Function declaration for modular code.
```c
#ifndef MATH_UTILS_H
#define MATH_UTILS_H

int solve_equation(int x);

#endif // MATH_UTILS_H
```

---

## 📄 `main/lib/math_utils.c`

Implements `solve_equation` using macro values to alter logic.
```c
#include "definitions.h"
#include "math_utils.h"

int solve_equation(int x) {
#if EQUATION_SYSTEM == 1
    return (USE_ADVANCED ? 3 * x + 5 : 2 * x + 1);
#elif EQUATION_SYSTEM == 2
    return (USE_ADVANCED ? x * x + 2 : x * x);
#else
    return 0;
#endif
}
```

---

## 📄 `main/main.c`

Entry point of the program that calls the function and prints the result.
```c
#include <stdio.h>
#include "lib/math_utils.h"

int main() {
    int input = 4;
    int result = solve_equation(input);
    printf("Result: %d\n", result);
    return 0;
}
```

---

## 🐍 `main/run.py`

Python script to automate:
- CMake configuration (with macro overrides)
- Compilation
- Execution
```python
import os
import subprocess

# --- Configuration ---
build_dir = os.path.join("main", "build")
macro1 = "EQUATION_SYSTEM"
macro1_val = 2

macro2 = "USE_ADVANCED"
macro2_val = 1

# --- Create build directory ---
os.makedirs(build_dir, exist_ok=True)

# --- Step 1: Configure with macros (override defaults if provided) ---
cmake_command = [
    "cmake",
    "-B", build_dir,
    "-S", "main",
    f"-DCMAKE_C_FLAGS=-D{macro1}={macro1_val} -D{macro2}={macro2_val}"
]

# --- Step 2: Build ---
build_command = ["cmake", "--build", build_dir]

# --- Step 3: Run the program ---
binary_path = os.path.join(build_dir, "my_program")

# --- Execute ---
print("Configuring...")
subprocess.run(cmake_command, check=True)

print("Building...")
subprocess.run(build_command, check=True)

print("Running program...")
subprocess.run([binary_path])
```

---

## ✅ How to Use

From the root folder, run:
```bash
python3 main/run.py
```

This will:
- Compile with `EQUATION_SYSTEM=2` and `USE_ADVANCED=1` (overriding the defaults in `definitions.h`)
- Run the resulting binary
- Print the result

You can modify the macro values in `run.py` or adjust the `CMake` logic as needed.

