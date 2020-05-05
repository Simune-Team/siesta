#!/bin/bash

# Run TS tests
MPI=${MPI:-mpirun -np 8}

# To run in serial mode, replace the 'mpirun' line
# by the appropriate incantation.

# Here we run the tests 
for d in 01-triclinic-aP.01-spg1 \
        01-triclinic-aP.02-spg2 \
        02-monoclinic-mP.01-spg3 \
        02-monoclinic-mP.02-spg4 \
        02-monoclinic-mS.03-spg5 \
        02-monoclinic-mP.04-spg6 \
        02-monoclinic-mP.05-spg7 \
        02-monoclinic-mS.06-spg8 \
        02-monoclinic-mS.07-spg9 \
        02-monoclinic-mP.08-spg10 \
        02-monoclinic-mP.09-spg11 \
        02-monoclinic-mS.10-spg12 \
        02-monoclinic-mP.11-spg13 \
        02-monoclinic-mP.12-spg14 \
        02-monoclinic-mS.13-spg15 \
        03-orthorhombic-oP.01-spg16 \
        03-orthorhombic-oP.02-spg17 \
        03-orthorhombic-oP.03-spg18 \
        03-orthorhombic-oP.04-spg19 \
        03-orthorhombic-oS.05-spg20 \
        03-orthorhombic-oS.06-spg21 \
        03-orthorhombic-oF.07-spg22 \
        03-orthorhombic-oI.08-spg23 \
        03-orthorhombic-oI.09-spg24 \
        03-orthorhombic-oP.10-spg25 \
        03-orthorhombic-oP.11-spg26 \
        03-orthorhombic-oP.12-spg27 \
        03-orthorhombic-oP.13-spg28 \
        03-orthorhombic-oP.14-spg29 \
        03-orthorhombic-oP.15-spg30 \
        03-orthorhombic-oP.16-spg31 \
        03-orthorhombic-oP.17-spg32 \
        03-orthorhombic-oP.18-spg33 \
        03-orthorhombic-oP.19-spg34 \
        03-orthorhombic-oS.20-spg35 \
        03-orthorhombic-oS.21-spg36 \
        03-orthorhombic-oS.22-spg37 \
        03-orthorhombic-oS.23-spg38 \
        03-orthorhombic-oS.24-spg39 \
        03-orthorhombic-oS.25-spg40 \
        03-orthorhombic-oS.26-spg41 \
        03-orthorhombic-oF.27-spg42 \
        03-orthorhombic-oF.28-spg43 \
        03-orthorhombic-oI.29-spg44 \
        03-orthorhombic-oI.30-spg45 \
        03-orthorhombic-oI.31-spg46 \
        03-orthorhombic-oP.32-spg47 \
        03-orthorhombic-oP.33-spg48 \
        03-orthorhombic-oP.34-spg49 \
        03-orthorhombic-oP.35-spg50 \
        03-orthorhombic-oP.36-spg51 \
        03-orthorhombic-oP.37-spg52 \
        03-orthorhombic-oP.38-spg53 \
        03-orthorhombic-oP.39-spg54 \
        03-orthorhombic-oP.40-spg55 \
        03-orthorhombic-oP.41-spg56 \
        03-orthorhombic-oP.42-spg57 \
        03-orthorhombic-oP.43-spg58 \
        03-orthorhombic-oP.44-spg59 \
        03-orthorhombic-oP.45-spg60 \
        03-orthorhombic-oP.46-spg61 \
        03-orthorhombic-oP.47-spg62 \
        03-orthorhombic-oS.48-spg63 \
        03-orthorhombic-oS.49-spg64 \
        03-orthorhombic-oS.50-spg65 \
        03-orthorhombic-oS.51-spg66 \
        03-orthorhombic-oS.52-spg67 \
        03-orthorhombic-oS.53-spg68 \
        03-orthorhombic-oF.54-spg69 \
        03-orthorhombic-oF.55-spg70 \
        03-orthorhombic-oI.56-spg71 \
        03-orthorhombic-oI.57-spg72 \
        03-orthorhombic-oI.58-spg73 \
        03-orthorhombic-oI.59-spg74 \
        04-tetragonal-oP.01-spg75 \
        04-tetragonal-oP.02-spg76 \
        04-tetragonal-oP.03-spg77 \
        04-tetragonal-oP.04-spg78 \
        04-tetragonal-tI.05-spg79 \
        04-tetragonal-tI.06-spg80 \
        04-tetragonal-oP.07-spg81 \
        04-tetragonal-tI.08-spg82 \
        04-tetragonal-oP.09-spg83 \
        04-tetragonal-oP.10-spg84 \
        04-tetragonal-oP.11-spg85 \
        04-tetragonal-oP.12-spg86 \
        04-tetragonal-tI.13-spg87 \
        04-tetragonal-tI.14-spg88 \
        04-tetragonal-oP.15-spg89 \
        04-tetragonal-oP.16-spg90 \
        04-tetragonal-oP.17-spg91 \
        04-tetragonal-oP.18-spg92 \
        04-tetragonal-oP.19-spg93 \
        04-tetragonal-oP.20-spg94 \
        04-tetragonal-oP.21-spg95 \
        04-tetragonal-oP.22-spg96 \
        04-tetragonal-tI.23-spg97 \
        04-tetragonal-tI.24-spg98 \
        04-tetragonal-oP.25-spg99 \
        04-tetragonal-oP.26-spg100 \
        04-tetragonal-oP.27-spg101 \
        04-tetragonal-oP.28-spg102 \
        04-tetragonal-oP.29-spg103 \
        04-tetragonal-oP.30-spg104 \
        04-tetragonal-oP.31-spg105 \
        04-tetragonal-oP.32-spg106 \
        04-tetragonal-tI.33-spg107 \
        04-tetragonal-tI.34-spg108 \
        04-tetragonal-tI.35-spg109 \
        04-tetragonal-tI.36-spg110 \
        04-tetragonal-oP.37-spg111 \
        04-tetragonal-oP.38-spg112 \
        04-tetragonal-oP.39-spg113 \
        04-tetragonal-oP.40-spg114 \
        04-tetragonal-oP.41-spg115 \
        04-tetragonal-oP.42-spg116 \
        04-tetragonal-oP.43-spg117 \
        04-tetragonal-oP.44-spg118 \
        04-tetragonal-tI.45-spg119 \
        04-tetragonal-tI.46-spg120 \
        04-tetragonal-tI.47-spg121 \
        04-tetragonal-tI.48-spg122 \
        04-tetragonal-oP.49-spg123 \
        04-tetragonal-oP.50-spg124 \
        04-tetragonal-oP.51-spg125 \
        04-tetragonal-oP.52-spg126 \
        04-tetragonal-oP.53-spg127 \
        04-tetragonal-oP.54-spg128 \
        04-tetragonal-oP.55-spg129 \
        04-tetragonal-oP.56-spg130 \
        04-tetragonal-oP.57-spg131 \
        04-tetragonal-oP.58-spg132 \
        04-tetragonal-oP.59-spg133 \
        04-tetragonal-oP.60-spg134 \
        04-tetragonal-oP.61-spg135 \
        04-tetragonal-oP.62-spg136 \
        04-tetragonal-oP.63-spg137 \
        04-tetragonal-oP.64-spg138 \
        04-tetragonal-tI.65-spg139 \
        04-tetragonal-tI.66-spg140 \
        04-tetragonal-tI.67-spg141 \
        04-tetragonal-tI.68-spg142 \
        05-hexagonal-hP.01-spg143 \
        05-hexagonal-hP.02-spg144 \
        05-hexagonal-hP.03-spg145 \
        05-hexagonal-hR.04-spg146 \
        05-hexagonal-hP.05-spg147 \
        05-hexagonal-hR.06-spg148 \
        05-hexagonal-hP.07-spg149 \
        05-hexagonal-hP.08-spg150 \
        05-hexagonal-hP.09-spg151 \
        05-hexagonal-hP.10-spg152 \
        05-hexagonal-hP.11-spg153 \
        05-hexagonal-hP.12-spg154 \
        05-hexagonal-hR.13-spg155 \
        05-hexagonal-hP.14-spg156 \
        05-hexagonal-hP.15-spg157 \
        05-hexagonal-hP.16-spg158 \
        05-hexagonal-hP.17-spg159 \
        05-hexagonal-hR.18-spg160 \
        05-hexagonal-hR.19-spg161 \
        05-hexagonal-hP.20-spg162 \
        05-hexagonal-hP.21-spg163 \
        05-hexagonal-hP.22-spg164 \
        05-hexagonal-hP.23-spg165 \
        05-hexagonal-hR.24-spg166 \
        05-hexagonal-hR.25-spg167 \
        05-hexagonal-hP.26-spg168 \
        05-hexagonal-hP.27-spg169 \
        05-hexagonal-hP.28-spg170 \
        05-hexagonal-hP.29-spg171 \
        05-hexagonal-hP.30-spg172 \
        05-hexagonal-hR.31-spg173 \
        05-hexagonal-hP.32-spg174 \
        05-hexagonal-hP.33-spg175 \
        05-hexagonal-hP.34-spg176 \
        05-hexagonal-hP.35-spg177 \
        05-hexagonal-hP.36-spg178 \
        05-hexagonal-hP.37-spg179 \
        05-hexagonal-hP.38-spg180 \
        05-hexagonal-hP.39-spg181 \
        05-hexagonal-hP.40-spg182 \
        05-hexagonal-hP.41-spg183 \
        05-hexagonal-hP.42-spg184 \
        05-hexagonal-hP.43-spg185 \
        05-hexagonal-hP.44-spg186 \
        05-hexagonal-hP.45-spg187 \
        05-hexagonal-hP.46-spg188 \
        05-hexagonal-hP.47-spg189 \
        05-hexagonal-hP.48-spg190 \
        05-hexagonal-hP.49-spg191 \
        05-hexagonal-hP.50-spg192 \
        05-hexagonal-hP.51-spg193 \
        05-hexagonal-hP.52-spg194 \
        06-cubic-cP.01-spg195 \
        06-cubic-cF.02-spg196 \
        06-cubic-cI.03-spg197 \
        06-cubic-cP.04-spg198 \
        06-cubic-cI.05-spg199 \
        06-cubic-cP.06-spg200 \
        06-cubic-cP.07-spg201 \
        06-cubic-cF.08-spg202 \
        06-cubic-cF.09-spg203 \
        06-cubic-cI.10-spg204 \
        06-cubic-cP.11-spg205 \
        06-cubic-cI.12-spg206 \
        06-cubic-cP.13-spg207 \
        06-cubic-cP.14-spg208 \
        06-cubic-cF.15-spg209 \
        06-cubic-cF.16-spg210 \
        06-cubic-cI.17-spg211 \
        06-cubic-cP.18-spg212 \
        06-cubic-cP.19-spg213 \
        06-cubic-cI.20-spg214 \
        06-cubic-cP.21-spg215 \
        06-cubic-cF.22-spg216 \
        06-cubic-cI.23-spg217 \
        06-cubic-cP.24-spg218 \
        06-cubic-cF.25-spg219 \
        06-cubic-cI.26-spg220 \
        06-cubic-cP.27-spg221 \
        06-cubic-cP.28-spg222 \
        06-cubic-cP.29-spg223 \
        06-cubic-cP.30-spg224 \
        06-cubic-cF.31-spg225 \
        06-cubic-cF.32-spg226 \
        06-cubic-cF.33-spg227 \
        06-cubic-cF.34-spg228 \
        06-cubic-cI.35-spg229 \
        06-cubic-cI.36-spg230 \

#        
       

do
    cd $d
    #make clean
    make MPI="$MPI"
    cd ..
done