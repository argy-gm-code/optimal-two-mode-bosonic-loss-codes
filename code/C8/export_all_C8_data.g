###############################################################################
# Export ALL irreps of the cyclic group C8.
###############################################################################
LoadPackage("RepnDecomp");;
LoadPackage("repsn");;

# --- 1) Pick the group -------------------------------------------------------
n := 8;;
groupTag := Concatenation("C", String(n));;
G := CyclicGroup(IsPermGroup, n);;

if Size(G) <> n then
    Error("Unexpected cyclic group order");
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

CharacterStrings := function(rep, elements)
    return List(
        elements,
        g -> String(TraceMat(Image(rep, g)))
    );
end;

FlatMatrixStrings := function(rep, elements)
    return List(
        elements,
        g -> ToFlatStrMat(Image(rep, g))
    );
end;

# --- 5) Export one JSON file per irrep --------------------------------------
for i in [1..Length(reps)] do
    r   := reps[i];
    dim := NrRows(r(One(G)));

    # Characters on all elements, in the ordering of 'elts'.
    char_row := CharacterStrings(r, elts);

    # Representation matrices for all elements, flattened row-wise.
    mats := FlatMatrixStrings(r, elts);

    # Matrices of the same fixed group generators in this irrep.
    genMats := FlatMatrixStrings(r, gens);

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
