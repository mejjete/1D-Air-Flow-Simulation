#!/usr/bin/awk -f

# Main AWK script
BEGIN {
    srand();  # Seed the random number generator
}

# Check if the first field is "head"
$1 == "\"head\"" {
    # Extract the numerical value
    gsub(/[^0-9]/, "", $3);
    # Generate a random number between 0 and 8
    random_val = int(rand() * 61);
    # Replace the numerical value with the random number
    $3 = random_val ",";
    $1 = "            " $1;
}

$1 == "\"tail\"" {
    # Extract the numerical value
    gsub(/[^0-9]/, "", $3);
    # Generate a random number between 0 and 8
    random_val = int(rand() * 61);
    # Replace the numerical value with the random number
    $3 = random_val ",";
    $1 = "            " $1 
}

# Print the modified or unmodified line
{
    print;
}