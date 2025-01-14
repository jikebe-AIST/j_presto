#!/bin/bash

# Default destination directory (user's bin directory)
DEFAULT_DEST_DIR="$HOME/local/bin"
# Default compiler
DEFAULT_COMPILER="ifort"

# Define required Python modules manually
REQUIRED_PYTHON_MODULES=("argparse" "os" "sys" "math" "re" "math" "subprocess" "shutil" "time" "itertools" "numpy" "statistics" "concurrent" "mpi4py")

# Check if required Python modules are installed
MISSING_MODULES=()
echo "Python module install check"
for module in "${REQUIRED_PYTHON_MODULES[@]}"; do
    if python3 -c "import $module" &> /dev/null; then
        echo "    $module is installed. OK"
    else
        MISSING_MODULES+=("$module")
    fi
done

if [ ${#MISSING_MODULES[@]} -ne 0 ]; then
    echo "The following Python modules are not installed: ${MISSING_MODULES[@]}"
    echo "Please install them using pip."
    exit 1
fi

# Prompt user for destination directory
echo
read -p "Destination directory (default: $DEFAULT_DEST_DIR): " user_input

# Use user input if provided, otherwise use default
if [ -n "$user_input" ]; then
    DEST_DIR="$user_input"
else
    DEST_DIR="$DEFAULT_DEST_DIR"
fi

# Prompt user for compiler
read -p "Fortran compiler (default: $DEFAULT_COMPILER): " user_compiler
echo

# Use user input if provided, otherwise use default
if [ -n "$user_compiler" ]; then
    COMPILER="$user_compiler"
else
    COMPILER="$DEFAULT_COMPILER"
fi

# Check if f2py is installed
if ! command -v f2py &> /dev/null; then
    echo "f2py is not installed."
    exit 1
fi

# Check if calc_Efield.f90 exists in the current directory
if [ ! -f "calc_Efield.f90" ]; then
    echo "calc_Efield.f90 does not exist in the current directory."
    exit 1
fi

# Execute the compile command
export CFLAGS="-std=c99"
f2py -c -m calc_Efield calc_Efield.f90 > /dev/null

# Check the result of the compilation
if [ $? -eq 0 ]; then
    echo "Module calc_Efield has been created successfully."

    # Get the generated file
    GENERATED_FILE=$(ls calc_Efield.cpython-*.so 2>/dev/null)

    # Check if the generated file exists
    if [ -z "$GENERATED_FILE" ]; then
        echo "The generated module file was not found."
        exit 1
    fi

    # Rename the file
    mv "$GENERATED_FILE" "calc_Efield.so"

    echo "The module has been renamed to calc_Efield.so."

    # Move the module to the existing sp directory
    SP_DIR="sp"

    if [ -d "$SP_DIR" ]; then
        mv "calc_Efield.so" "$SP_DIR/"
        echo "The module has been moved to the $SP_DIR directory."
    else
        echo "The sp directory does not exist in the current directory."
        exit 1
    fi

    # Create the J_PRESTO directory
    PRESTO_DIR="J_PRESTO"
    mkdir -p "$PRESTO_DIR"

    # Copy the sp directory to the J_PRESTO directory
    cp -r "$SP_DIR" "$PRESTO_DIR"
    
    echo "The $SP_DIR directory has been copied to $PRESTO_DIR."

    # Copy other necessary files to the J_PRESTO directory
    cp LICENSE.md j_presto* manual.txt "$PRESTO_DIR"
    cp -r template "$PRESTO_DIR"
    cp -r database "$PRESTO_DIR"

    echo "Other files have been copied to $PRESTO_DIR."

    # Store the current directory
    CURRENT_DIR=$(pwd)

    build_and_copy() {
        local src_dir=$1
	local exe_names=($2)

        if [ -d "$src_dir" ]; then
            cd "$src_dir" || exit

            # Replace the FC variable in the Makefile with the chosen compiler
            sed -i "s/^FC = .*/FC = $COMPILER/" Makefile

            if make; then
                echo "Make completed successfully in $src_dir."
                # Copy the executable to the J_PRESTO directory
		for exe_name in "${exe_names[@]}"; do
                    cp "$exe_name" "../$PRESTO_DIR/sp/"
                    echo "Executable $exe_name has been copied to $PRESTO_DIR."
                done
            else
                echo "Make failed in $src_dir."
                exit 1
            fi

            # Return to the original directory
            cd "$CURRENT_DIR" || exit
        else
            echo "$src_dir directory does not exist."
            exit 1
        fi
    }

    # Build and copy executables
    build_and_copy "src_md_run" "md_run mkshkl"
    build_and_copy "src_GEprep" "GEprep"
    build_and_copy "src_Ens_Ana" "Ens_Ana"

    # Move the J_PRESTO directory to the specified destination
    # Check if the destination already has a J_PRESTO directory
    if [ -d "$DEST_DIR/J_PRESTO" ]; then
        echo "Old J_PRESTO directory found in $DEST_DIR. Deleting it."
        rm -rf "$DEST_DIR/J_PRESTO"
    fi

    # Move the new J_PRESTO directory to the destination
    mv "$PRESTO_DIR" "$DEST_DIR" && echo "The J_PRESTO directory has been moved to $DEST_DIR." || echo "Move failed."

    # Output environment variable setup messages
    echo -e "\n\n"  # Add empty lines for spacing
    echo "!! Please add the following lines to your .bashrc or .bash_profile to set up the environment variables:"
    echo "    export J_PRESTO_PATH=\"$DEST_DIR/$PRESTO_DIR\""
    echo "    export J_PRESTO_AMBER_DATABASE_PATH=\"<path_to_your_amber_database>\""
    echo "    export PATH=\"\$PATH:\$J_PRESTO_PATH\""
    echo
else
    echo "An error occurred during the creation of the module."
    exit 1
fi

