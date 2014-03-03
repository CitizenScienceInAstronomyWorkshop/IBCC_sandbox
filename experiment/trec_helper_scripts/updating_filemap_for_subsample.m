To select a subsample of the TREC data.

1. Select the file map indexes to keep from the original filemap. Remove the other rows
2. Create new file map indexes and store a map from old to new.
3. Change the indexes in the filemap to the new set.
4. Select the rows from the feature file that correspond to the selected indices.
5. Switch the selected indexes in the feature file to new file map indexes.
6. Adjust the size parameters at start of feature file. 

Alternatively:

Use the full set of document IDs, the original filemap as for TREC.

1. Select the file map indexes to keep from the original filemap.
2. Write the selected indexes to file.
3. Select the rows from the feature file that correspond to the selected indices. 
a) Remove the others. b) Keep the old filemap indexes
4. Add in a compression code before running the classifier: removes empty lines in the feature matrix.
5. Add in a line of code to map the selected documents and results back to original filemap indexes.

