import React, { useState } from 'react';
import {
  Box,
  Button,
  Card,
  CardContent,
  CardMedia,
  Container,
  Grid,
  Paper,
  TextField,
  Typography,
  Slider,
  Alert,
  CircularProgress,
  Chip,
  Divider,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
} from '@mui/material';
import {
  Search as SearchIcon,
  Science as ScienceIcon,
} from '@mui/icons-material';
import drugDiscoveryApi from '../services/api';

function SearchPage() {
  // State for form inputs
  const [smiles, setSmiles] = useState('');
  const [numResults, setNumResults] = useState(10);
  
  // State for API response
  const [similarMolecules, setSimilarMolecules] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [success, setSuccess] = useState(false);
  
  // Handle form submission
  const handleSubmit = async (e) => {
    e.preventDefault();
    setLoading(true);
    setError(null);
    setSuccess(false);
    
    try {
      const response = await drugDiscoveryApi.searchSimilarMolecules(
        smiles,
        numResults
      );
      
      setSimilarMolecules(response.similar_molecules);
      setSuccess(true);
      window.scrollTo({ top: document.getElementById('results').offsetTop, behavior: 'smooth' });
    } catch (error) {
      setError('Error searching for similar molecules. Please try again.');
      console.error(error);
    } finally {
      setLoading(false);
    }
  };
  
  // Example molecules for quick selection
  const exampleMolecules = [
    { name: 'Aspirin', smiles: 'CC(=O)Oc1ccccc1C(=O)O' },
    { name: 'Ibuprofen', smiles: 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O' },
    { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' },
    { name: 'Paracetamol', smiles: 'CC(=O)NC1=CC=C(C=C1)O' },
  ];
  
  return (
    <Container maxWidth="lg">
      <Typography variant="h4" component="h1" gutterBottom sx={{ mb: 3 }}>
        <SearchIcon sx={{ mr: 1, verticalAlign: 'middle' }} />
        Search Similar Molecules
      </Typography>
      
      <Typography variant="body1" paragraph>
        Enter a SMILES string to search for similar molecules in our database. The search uses molecular fingerprints to find structurally similar compounds.
      </Typography>
      
      <Grid container spacing={4}>
        {/* Form Section */}
        <Grid item xs={12} md={5}>
          <Paper elevation={3} sx={{ p: 3, borderRadius: 2 }}>
            <form onSubmit={handleSubmit}>
              <Typography variant="h6" gutterBottom>
                Search Parameters
              </Typography>
              
              <Box sx={{ mb: 3 }}>
                <TextField
                  label="SMILES String"
                  value={smiles}
                  onChange={(e) => setSmiles(e.target.value)}
                  fullWidth
                  required
                  placeholder="Enter SMILES notation"
                  helperText="e.g., CC(=O)Oc1ccccc1C(=O)O for Aspirin"
                />
              </Box>
              
              <Typography gutterBottom>
                Example Molecules:
              </Typography>
              <Box sx={{ mb: 3, display: 'flex', flexWrap: 'wrap', gap: 1 }}>
                {exampleMolecules.map((mol, index) => (
                  <Chip
                    key={index}
                    label={mol.name}
                    onClick={() => setSmiles(mol.smiles)}
                    clickable
                  />
                ))}
              </Box>
              
              <Box sx={{ mb: 3 }}>
                <Typography gutterBottom>
                  Number of Results
                </Typography>
                <Slider
                  value={numResults}
                  onChange={(e, value) => setNumResults(value)}
                  valueLabelDisplay="auto"
                  step={1}
                  marks
                  min={1}
                  max={20}
                />
              </Box>
              
              <Button
                type="submit"
                variant="contained"
                color="primary"
                fullWidth
                size="large"
                disabled={loading || !smiles}
                startIcon={loading ? <CircularProgress size={24} /> : <SearchIcon />}
              >
                {loading ? 'Searching...' : 'Search Similar Molecules'}
              </Button>
              
              {error && (
                <Alert severity="error" sx={{ mt: 2 }}>
                  {error}
                </Alert>
              )}
              
              {success && (
                <Alert severity="success" sx={{ mt: 2 }}>
                  Found {similarMolecules.length} similar molecules!
                </Alert>
              )}
            </form>
          </Paper>
        </Grid>
        
        {/* Results Section */}
        <Grid item xs={12} md={7} id="results">
          <Paper elevation={3} sx={{ p: 3, borderRadius: 2 }}>
            <Typography variant="h6" gutterBottom>
              Similar Molecules
            </Typography>
            
            {similarMolecules.length === 0 ? (
              <Box sx={{ textAlign: 'center', py: 4 }}>
                <SearchIcon sx={{ fontSize: 60, color: 'text.secondary', opacity: 0.3 }} />
                <Typography variant="body1" color="text.secondary" sx={{ mt: 2 }}>
                  No similar molecules found yet. Enter a SMILES string and click "Search Similar Molecules" to start.
                </Typography>
              </Box>
            ) : (
              <TableContainer>
                <Table>
                  <TableHead>
                    <TableRow>
                      <TableCell>Name</TableCell>
                      <TableCell>Structure</TableCell>
                      <TableCell>SMILES</TableCell>
                      <TableCell>Similarity</TableCell>
                      <TableCell>Target</TableCell>
                    </TableRow>
                  </TableHead>
                  <TableBody>
                    {similarMolecules.map((molecule, index) => (
                      <TableRow key={index}>
                        <TableCell>{molecule.name}</TableCell>
                        <TableCell>
                          <Box
                            component="img"
                            sx={{
                              height: 100,
                              width: 100,
                              objectFit: 'contain',
                            }}
                            alt={molecule.name}
                            src={`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(molecule.smiles)}/PNG`}
                          />
                        </TableCell>
                        <TableCell>
                          <Typography variant="body2" sx={{ fontFamily: 'monospace', fontSize: '0.8rem' }}>
                            {molecule.smiles.length > 20 ? `${molecule.smiles.substring(0, 20)}...` : molecule.smiles}
                          </Typography>
                        </TableCell>
                        <TableCell>
                          <Chip 
                            label={`${(molecule.similarity * 100).toFixed(1)}%`} 
                            color={molecule.similarity > 0.7 ? 'success' : molecule.similarity > 0.5 ? 'primary' : 'default'}
                          />
                        </TableCell>
                        <TableCell>{molecule.target}</TableCell>
                      </TableRow>
                    ))}
                  </TableBody>
                </Table>
              </TableContainer>
            )}
          </Paper>
        </Grid>
      </Grid>
    </Container>
  );
}

export default SearchPage; 