###############################################################################
# Export ALL irreps of SL(2,5) = 2I in a dimension-agnostic format.
###############################################################################
LoadPackage("RepnDecomp");;
LoadPackage("repsn");;

# --- 1) Pick your group ------------------------------------------------------
G := SL(2,5);; # For 2I 

#G := SmallGroup(48, 28);;; #for 2O 
#G := SmallGroup(24,3);;  # for 2T 

# --- 2) All irreps -----------------------------------------------------------
reps := IrreducibleRepresentations(G);;
dims := List(reps, r -> NrRows(r(One(G))));
Print("All irrep dimensions: ", dims, "\n");
Print("Total irreps: ", Length(reps), "\n");

# --- 3) Conjugacy classes and elements --------------------------------------
classes    := ConjugacyClasses(G);;
classSizes := List(classes, Size);;
elts       := Elements(G);;
gens       := GeneratorsOfGroup(G);;
genPositions := List(gens, g -> Position(elts, g));;

Print("Generator positions in element list: ", genPositions, "\n");

# --- 4) Helper functions -----------------------------------------------------
ToFlatStrMat := function(M)
    # flatten d×d matrix row-wise to list of strings
    return List(Flat(M), x -> String(x));
end;

WriteText := function(filename, txt)
    local f;
    f := OutputTextFile(filename, false);;
    SetPrintFormattingStatus(f, false);
    PrintTo(f, txt, "\n");
    CloseStream(f);
end;

# --- 5) Export one file per irrep -------------------------------------------
for i in [1..Length(reps)] do
    r   := reps[i];
    dim := NrRows(r(One(G)));

    # Characters on all elements (same order as 'elts')
    char_row := List(elts, g -> String(TraceMat(Image(r, g))));

    # Matrices for all elements (flattened row-wise)
    mats := List(elts, g -> ToFlatStrMat(Image(r, g)));

    # Matrices for the fixed group generators above, in this same irrep.
    # These are not recomputed from a separate representation.
    genMats := List(gens, g -> ToFlatStrMat(Image(r, g)));

    # File name encodes dimension and index, same pattern as before:
    # e.g. "irrep_dim2_1.json", "irrep_dim3_2.json", etc.
    filename := Concatenation("irrep_dim", String(dim), "_2I_", String(i), ".json");
    Print("Writing ", filename, " ...\n");

    txt := "{";
    Append(txt, "\"dimension\":");    Append(txt, String(dim));        Append(txt, ",");
    Append(txt, "\"class_sizes\":");  Append(txt, String(classSizes)); Append(txt, ",");
    Append(txt, "\"characters\":");   Append(txt, String(char_row));   Append(txt, ",");
    Append(txt, "\"matrices\":");     Append(txt, String(mats));       Append(txt, ",");
    Append(txt, "\"generator_indices\":");  Append(txt, String(genPositions)); Append(txt, ",");
    Append(txt, "\"generator_matrices\":"); Append(txt, String(genMats));
    Append(txt, "}");

    WriteText(filename, txt);
od;

Print("✓ Done. One JSON file per irrep of the selected group.\n");
###############################################################################
