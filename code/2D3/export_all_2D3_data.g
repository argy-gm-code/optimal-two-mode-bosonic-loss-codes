###############################################################################
# Export ALL irreps of the binary dihedral group 2D3 = Dic_3.
###############################################################################
LoadPackage("RepnDecomp");;
LoadPackage("repsn");;

# --- 1) Pick the group -------------------------------------------------------
n := 3;;
groupTag := Concatenation("2D", String(n));;

# Binary dihedral group of order 4n in its natural two-dimensional
# realization inside SU(2).
G := DicyclicGroup(IsMatrixGroup, CF(2*n), 4*n);;

if Size(G) <> 4*n then
    Error("Unexpected binary-dihedral group order");
fi;

if not ForAll(GeneratorsOfGroup(G), g -> DeterminantMat(g) = 1) then
    Error("The natural group generators are not in SU(2)");
fi;

# --- 2) All irreps -----------------------------------------------------------
reps := IrreducibleRepresentations(G);;
dims := List(reps, r -> NrRows(r(One(G))));;

Print("Group: ", groupTag, "\n");
Print("Group order: ", Size(G), "\n");
Print("All irrep dimensions: ", dims, "\n");
Print("Total irreps: ", Length(reps), "\n");

if Sum(dims, d -> d^2) <> Size(G) then
    Error("Irrep dimensions do not satisfy the sum-of-squares identity");
fi;

# --- 3) Conjugacy classes and elements --------------------------------------
classes      := ConjugacyClasses(G);;
classSizes   := List(classes, Size);;
elts         := Elements(G);;
gens         := GeneratorsOfGroup(G);;
genPositions := List(gens, g -> Position(elts, g));;

Print("Generator positions in element list: ", genPositions, "\n");

# Identify the exported irrep equivalent to the natural SU(2) representation.
naturalCharacter := List(elts, g -> TraceMat(g));;
fundamentalIndices := Filtered(
    [1..Length(reps)],
    i -> dims[i] = 2 and
         ForAll(
             [1..Length(elts)],
             j -> TraceMat(Image(reps[i], elts[j])) = naturalCharacter[j]
         )
);;

if Length(fundamentalIndices) <> 1 then
    Error("Could not uniquely identify the natural SU(2) irrep");
fi;

fundamentalIndex := fundamentalIndices[1];;
Print(
    "Physical fundamental irrep file: irrep_dim2_",
    groupTag, "_", fundamentalIndex, ".json\n"
);

# --- 4) Helper functions -----------------------------------------------------
ToFlatStrMat := function(M)
    # Flatten a d-by-d matrix row-wise to a list of strings.
    return List(Flat(M), x -> String(x));
end;

WriteText := function(filename, txt)
    local f;
    f := OutputTextFile(filename, false);;
    SetPrintFormattingStatus(f, false);
    PrintTo(f, txt, "\n");
    CloseStream(f);
end;

# --- 5) Export one JSON file per irrep --------------------------------------
for i in [1..Length(reps)] do
    r   := reps[i];
    dim := NrRows(r(One(G)));

    # Characters on all elements, in the ordering of 'elts'.
    char_row := List(elts, g -> String(TraceMat(Image(r, g))));

    # Representation matrices for all elements, flattened row-wise.
    mats := List(elts, g -> ToFlatStrMat(Image(r, g)));

    # Matrices of the same fixed group generators in this irrep.
    genMats := List(gens, g -> ToFlatStrMat(Image(r, g)));

    filename := Concatenation(
        "irrep_dim", String(dim), "_", groupTag, "_", String(i), ".json"
    );
    Print("Writing ", filename, " ...\n");

    txt := "{";
    Append(txt, "\"dimension\":");
    Append(txt, String(dim));
    Append(txt, ",");
    Append(txt, "\"class_sizes\":");
    Append(txt, String(classSizes));
    Append(txt, ",");
    Append(txt, "\"characters\":");
    Append(txt, String(char_row));
    Append(txt, ",");
    Append(txt, "\"matrices\":");
    Append(txt, String(mats));
    Append(txt, ",");
    Append(txt, "\"generator_indices\":");
    Append(txt, String(genPositions));
    Append(txt, ",");
    Append(txt, "\"generator_matrices\":");
    Append(txt, String(genMats));
    Append(txt, "}");

    WriteText(filename, txt);
od;

Print("Done. One JSON file per irrep of ", groupTag, ".\n");
###############################################################################
