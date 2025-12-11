#!/bin/bash

#!/bin/bash

# Verificación de argumento
if [ $# -lt 1 ]; then
    echo "Uso: $0 <ejecutable>"
    exit 1
fi

EXEC=$1

# Arrays solicitados
PROCS=(2 4 8 16 24 32 40)
SIZES=(100 500 750 1000)

# Número de repeticiones por caso
REPS=5

for np in "${PROCS[@]}"; do
    echo "==============================="
    echo " Ejecutando con np = $np"
    echo "==============================="

    for size in "${SIZES[@]}"; do
        echo "---- Tamaño: $size x $size ----"

        for ((i = 1; i <= REPS; i++)); do
            echo "Ejecutando iteración $i..."
            mpirun -np "$np" "$EXEC" "$size" "$size"
        done

        echo ""
    done
done


