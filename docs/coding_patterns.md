# NeoPZ Coding Patterns and Best Practices

## TPZAdmChunkVector: Reference vs Copy

When working with `TPZAdmChunkVector` (especially when accessing element vectors via `ElementVec()`), it's important to understand when to use a reference versus when to make a copy.

### Reference vs Copy Syntax

```cpp
// REFERENCE (using &) - No copy is made, just an alias to the original vector
TPZAdmChunkVector<TPZGeoEl *> &elvec = cmesh.Reference()->ElementVec();

// COPY (without &) - A complete copy of the vector is created
TPZAdmChunkVector<TPZGeoEl *> elvec = cmesh.Reference()->ElementVec();
```

### When to Use References (`&`)

Use a **reference** when:
- You only need to **read** from the vector
- You iterate through elements without modifying the vector structure
- Performance is critical and you want to avoid copying overhead

**Example locations:**
- `Mesh/pzcreateapproxspace.cpp:253` - BuildMesh() only reads geometric elements
- `Mesh/pzcmesh.cpp:420` - AutoBuildContDisc() only reads elements
- `Pre/pzbuildmultiphysicsmesh.cpp:452` - BuildHybridMesh() iterates for reading

```cpp
void TPZCreateApproximationSpace::BuildMesh(TPZCompMesh &cmesh, const std::set<int> &MaterialIDs) const {
    // Safe to use reference - only reading from the vector
    TPZAdmChunkVector<TPZGeoEl *> &elvec = cmesh.Reference()->ElementVec();
    for(int64_t i=0; i<elvec.NElements(); i++) {
        TPZGeoEl *gel = elvec[i];
        // Only reading, not modifying the vector
    }
}
```

### When to Use Copies (without `&`)

Use a **copy** when:
- You need to **delete elements** from the mesh during iteration
- Operations might **add new elements** to the mesh (like refinement/division)
- You need to iterate over the **original state** while the mesh structure changes

**Example locations:**
- `Pre/pzbuildmultiphysicsmesh.cpp:660` - Deletes interface elements during iteration
- `Pre/pzbuildmultiphysicsmesh.cpp:677` - Divides elements which adds new elements to mesh

```cpp
void TPZBuildMultiphysicsMesh::UniformRefineCompMesh(TPZCompMesh *cMesh, int ndiv, bool isLagrMult) {
    // MUST use a copy - elements are deleted during iteration
    TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
    for(int64_t el=0; el < elvec.NElements(); el++){
        TPZCompEl * compEl = elvec[el];
        if(compEl && compEl->IsInterface()){
            delete compEl;  // Modifies the original mesh
        }
    }
}
```

### Common Misconceptions

❌ **INCORRECT**: "The `&` creates a copy"
- The `&` symbol means **reference**, not copy
- It creates an alias to the original vector
- No memory duplication occurs

✅ **CORRECT**: "The `&` creates a reference"
- Changes to the referenced vector affect the original
- No copying overhead
- Must be careful if vector structure changes during iteration

### Performance Considerations

- **References**: Zero copy overhead, very efficient
- **Copies**: Full vector copy, includes all pointers but creates a snapshot
- Rule of thumb: Use references unless you have a specific reason to copy

### Safety Guidelines

1. **For read-only operations**: Always use references (`&`)
2. **For element deletion**: Use a copy to avoid iterator invalidation
3. **For mesh refinement/division**: Use a copy since new elements are added
4. **For element creation**: Usually safe with reference if you don't iterate after creation

