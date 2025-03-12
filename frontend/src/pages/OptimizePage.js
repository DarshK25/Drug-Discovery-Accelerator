import React, { useState } from 'react';
import {
  Box,
  Button,
  Card,
  CardContent,
  CardMedia,
  Container,
  Divider,
  FormControl,
  FormControlLabel,
  FormGroup,
  Grid,
  InputLabel,
  MenuItem,
  Paper,
  Select,
  Slider,
  Switch,
  TextField,
  Typography,
  Alert,
  CircularProgress,
  Chip,
  Tooltip,
  IconButton,
  Tabs,
  Tab,
} from '@mui/material';
import {
  Tune as TuneIcon,
  Science as ScienceIcon,
  Compare as CompareIcon,
  Download as DownloadIcon,
  Bookmark as BookmarkIcon,
} from '@mui/icons-material';
import drugDiscoveryApi from '../services/api';

function OptimizePage() {
  // State for form inputs
  const [smiles, setSmiles] = useState('');
  const [targetProperties, setTargetProperties] = useState({
    logP: { target: 2.5, weight: 1.0 },
    solubility: { target: 'High', weight: 1.0 },
    qed: { target: 0.8, weight: 1.0 },
    synthetic_accessibility: { target: 3.0, weight: 1.0 },
  });
  
  // State for API response
  const [optimizedMolecules, setOptimizedMolecules] = useState([]);
  const [originalProperties, setOriginalProperties] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [success, setSuccess] = useState(false);
  
  // State for tabs
  const [tabValue, setTabValue] = useState(0);
  
  // Handle form submission
  const handleSubmit = async (e) => {
    e.preventDefault();
    setLoading(true);
    setError(null);
    setSuccess(false);
    
    try {
      // First, get properties of the original molecule
      const propertiesResponse = await drugDiscoveryApi.predictProperties(smiles);
      setOriginalProperties(propertiesResponse.properties);
      
      // Then, optimize the molecule
      const response = await drugDiscoveryApi.optimizeMolecule(
        smiles,
        targetProperties
      );
      
      setOptimizedMolecules(response.optimized_molecules);
      setSuccess(true);
      window.scrollTo({ top: document.getElementById('results').offsetTop, behavior: 'smooth' });
    } catch (error) {
      setError('Error optimizing molecule. Please try again.');
      console.error(error);
    } finally {
      setLoading(false);
    }
  };
  
  // Handle property changes
  const handlePropertyChange = (property, type, value) => {
    setTargetProperties((prev) => ({
      ...prev,
      [property]: {
        ...prev[property],
        [type]: value,
      },
    }));
  };
  
  // Handle tab change
  const handleTabChange = (event, newValue) => {
    setTabValue(newValue);
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
        <TuneIcon sx={{ mr: 1, verticalAlign: 'middle' }} />
        Optimize Molecules
      </Typography>
      
      <Typography variant="body1" paragraph>
        Enter a SMILES string and specify target properties to optimize the molecule. Our AI will generate variations that better match your desired properties.
      </Typography>
      
      <Grid container spacing={4}>
        {/* Form Section */}
        <Grid item xs={12} md={5}>
          <Paper elevation={3} sx={{ p: 3, borderRadius: 2 }}>
            <form onSubmit={handleSubmit}>
              <Typography variant="h6" gutterBottom>
                Starting Molecule
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
              
              <Divider sx={{ my: 3 }} />
              
              <Typography variant="h6" gutterBottom>
                Target Properties
              </Typography>
              
              <Box sx={{ mb: 3 }}>
                <Typography gutterBottom>
                  LogP (Lipophilicity) - Target: {targetProperties.logP.target}
                </Typography>
                <Slider
                  value={targetProperties.logP.target}
                  onChange={(e, value) => handlePropertyChange('logP', 'target', value)}
                  valueLabelDisplay="auto"
                  step={0.1}
                  marks={[
                    { value: 0, label: '0' },
                    { value: 2.5, label: '2.5' },
                    { value: 5, label: '5' },
                  ]}
                  min={0}
                  max={5}
                />
                <Typography gutterBottom sx={{ mt: 2 }}>
                  Importance Weight: {targetProperties.logP.weight}
                </Typography>
                <Slider
                  value={targetProperties.logP.weight}
                  onChange={(e, value) => handlePropertyChange('logP', 'weight', value)}
                  valueLabelDisplay="auto"
                  step={0.1}
                  min={0}
                  max={2}
                />
              </Box>
              
              <Box sx={{ mb: 3 }}>
                <Typography gutterBottom>
                  QED (Drug-likeness) - Target: {targetProperties.qed.target}
                </Typography>
                <Slider
                  value={targetProperties.qed.target}
                  onChange={(e, value) => handlePropertyChange('qed', 'target', value)}
                  valueLabelDisplay="auto"
                  step={0.05}
                  marks={[
                    { value: 0, label: '0' },
                    { value: 0.5, label: '0.5' },
                    { value: 1, label: '1' },
                  ]}
                  min={0}
                  max={1}
                />
                <Typography gutterBottom sx={{ mt: 2 }}>
                  Importance Weight: {targetProperties.qed.weight}
                </Typography>
                <Slider
                  value={targetProperties.qed.weight}
                  onChange={(e, value) => handlePropertyChange('qed', 'weight', value)}
                  valueLabelDisplay="auto"
                  step={0.1}
                  min={0}
                  max={2}
                />
              </Box>
              
              <Box sx={{ mb: 3 }}>
                <Typography gutterBottom>
                  Synthetic Accessibility - Target: {targetProperties.synthetic_accessibility.target}
                </Typography>
                <Slider
                  value={targetProperties.synthetic_accessibility.target}
                  onChange={(e, value) => handlePropertyChange('synthetic_accessibility', 'target', value)}
                  valueLabelDisplay="auto"
                  step={0.5}
                  marks={[
                    { value: 1, label: 'Easy (1)' },
                    { value: 5, label: 'Medium (5)' },
                    { value: 10, label: 'Hard (10)' },
                  ]}
                  min={1}
                  max={10}
                />
                <Typography gutterBottom sx={{ mt: 2 }}>
                  Importance Weight: {targetProperties.synthetic_accessibility.weight}
                </Typography>
                <Slider
                  value={targetProperties.synthetic_accessibility.weight}
                  onChange={(e, value) => handlePropertyChange('synthetic_accessibility', 'weight', value)}
                  valueLabelDisplay="auto"
                  step={0.1}
                  min={0}
                  max={2}
                />
              </Box>
              
              <Box sx={{ mb: 3 }}>
                <FormControl fullWidth>
                  <InputLabel>Target Solubility</InputLabel>
                  <Select
                    value={targetProperties.solubility.target}
                    onChange={(e) => handlePropertyChange('solubility', 'target', e.target.value)}
                    label="Target Solubility"
                  >
                    <MenuItem value="Low">Low</MenuItem>
                    <MenuItem value="Moderate">Moderate</MenuItem>
                    <MenuItem value="High">High</MenuItem>
                  </Select>
                </FormControl>
                <Typography gutterBottom sx={{ mt: 2 }}>
                  Importance Weight: {targetProperties.solubility.weight}
                </Typography>
                <Slider
                  value={targetProperties.solubility.weight}
                  onChange={(e, value) => handlePropertyChange('solubility', 'weight', value)}
                  valueLabelDisplay="auto"
                  step={0.1}
                  min={0}
                  max={2}
                />
              </Box>
              
              <Button
                type="submit"
                variant="contained"
                color="primary"
                fullWidth
                size="large"
                disabled={loading || !smiles}
                startIcon={loading ? <CircularProgress size={24} /> : <TuneIcon />}
              >
                {loading ? 'Optimizing...' : 'Optimize Molecule'}
              </Button>
              
              {error && (
                <Alert severity="error" sx={{ mt: 2 }}>
                  {error}
                </Alert>
              )}
              
              {success && (
                <Alert severity="success" sx={{ mt: 2 }}>
                  Successfully optimized molecule with {optimizedMolecules.length} variations!
                </Alert>
              )}
            </form>
          </Paper>
        </Grid>
        
        {/* Results Section */}
        <Grid item xs={12} md={7} id="results">
          <Paper elevation={3} sx={{ p: 3, borderRadius: 2 }}>
            <Typography variant="h6" gutterBottom>
              Optimization Results
            </Typography>
            
            {optimizedMolecules.length === 0 ? (
              <Box sx={{ textAlign: 'center', py: 4 }}>
                <TuneIcon sx={{ fontSize: 60, color: 'text.secondary', opacity: 0.3 }} />
                <Typography variant="body1" color="text.secondary" sx={{ mt: 2 }}>
                  No optimized molecules yet. Enter a SMILES string and target properties, then click "Optimize Molecule" to start.
                </Typography>
              </Box>
            ) : (
              <>
                <Box sx={{ borderBottom: 1, borderColor: 'divider', mb: 2 }}>
                  <Tabs value={tabValue} onChange={handleTabChange} aria-label="molecule tabs">
                    <Tab label="Original" />
                    {optimizedMolecules.map((_, index) => (
                      <Tab key={index} label={`Variant ${index + 1}`} />
                    ))}
                  </Tabs>
                </Box>
                
                {/* Original Molecule */}
                {tabValue === 0 && originalProperties && (
                  <Card>
                    <Grid container>
                      <Grid item xs={12} sm={4}>
                        <CardMedia
                          component="img"
                          sx={{ height: 200 }}
                          image={`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(smiles)}/PNG`}
                          alt="Original Molecule"
                        />
                      </Grid>
                      <Grid item xs={12} sm={8}>
                        <CardContent>
                          <Typography variant="h6" gutterBottom>
                            Original Molecule
                            <Tooltip title="Save molecule">
                              <IconButton size="small" sx={{ ml: 1 }}>
                                <BookmarkIcon fontSize="small" />
                              </IconButton>
                            </Tooltip>
                            <Tooltip title="Download structure">
                              <IconButton size="small">
                                <DownloadIcon fontSize="small" />
                              </IconButton>
                            </Tooltip>
                          </Typography>
                          
                          <Typography variant="body2" color="text.secondary" sx={{ mb: 2, fontFamily: 'monospace' }}>
                            {smiles}
                          </Typography>
                          
                          <Grid container spacing={1}>
                            <Grid item xs={6}>
                              <Typography variant="body2" color="text.secondary">
                                LogP:
                              </Typography>
                              <Typography variant="body2">
                                {originalProperties.logP}
                              </Typography>
                            </Grid>
                            <Grid item xs={6}>
                              <Typography variant="body2" color="text.secondary">
                                QED:
                              </Typography>
                              <Typography variant="body2">
                                {originalProperties.qed}
                              </Typography>
                            </Grid>
                            <Grid item xs={6}>
                              <Typography variant="body2" color="text.secondary">
                                Synthetic Accessibility:
                              </Typography>
                              <Typography variant="body2">
                                {originalProperties.synthetic_accessibility}
                              </Typography>
                            </Grid>
                            <Grid item xs={6}>
                              <Typography variant="body2" color="text.secondary">
                                Solubility:
                              </Typography>
                              <Typography variant="body2">
                                {originalProperties.solubility}
                              </Typography>
                            </Grid>
                          </Grid>
                        </CardContent>
                      </Grid>
                    </Grid>
                  </Card>
                )}
                
                {/* Optimized Molecules */}
                {optimizedMolecules.map((molecule, index) => (
                  tabValue === index + 1 && (
                    <Card key={index}>
                      <Grid container>
                        <Grid item xs={12} sm={4}>
                          <CardMedia
                            component="img"
                            sx={{ height: 200 }}
                            image={molecule.image_url || `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(molecule.smiles)}/PNG`}
                            alt={`Optimized Molecule ${index + 1}`}
                          />
                        </Grid>
                        <Grid item xs={12} sm={8}>
                          <CardContent>
                            <Typography variant="h6" gutterBottom>
                              Optimized Variant {index + 1}
                              <Tooltip title="Save molecule">
                                <IconButton size="small" sx={{ ml: 1 }}>
                                  <BookmarkIcon fontSize="small" />
                                </IconButton>
                              </Tooltip>
                              <Tooltip title="Download structure">
                                <IconButton size="small">
                                  <DownloadIcon fontSize="small" />
                                </IconButton>
                              </Tooltip>
                              <Tooltip title="Compare with original">
                                <IconButton size="small">
                                  <CompareIcon fontSize="small" />
                                </IconButton>
                              </Tooltip>
                            </Typography>
                            
                            <Typography variant="body2" color="text.secondary" sx={{ mb: 2, fontFamily: 'monospace' }}>
                              {molecule.smiles}
                            </Typography>
                            
                            <Grid container spacing={1}>
                              <Grid item xs={6}>
                                <Typography variant="body2" color="text.secondary">
                                  LogP:
                                </Typography>
                                <Typography variant="body2">
                                  {molecule.properties.logP}
                                  {originalProperties && (
                                    <Chip 
                                      size="small" 
                                      label={`${(molecule.properties.logP - originalProperties.logP).toFixed(2)}`}
                                      color={Math.abs(molecule.properties.logP - targetProperties.logP.target) < Math.abs(originalProperties.logP - targetProperties.logP.target) ? 'success' : 'error'}
                                      sx={{ ml: 1 }}
                                    />
                                  )}
                                </Typography>
                              </Grid>
                              <Grid item xs={6}>
                                <Typography variant="body2" color="text.secondary">
                                  QED:
                                </Typography>
                                <Typography variant="body2">
                                  {molecule.properties.qed}
                                  {originalProperties && (
                                    <Chip 
                                      size="small" 
                                      label={`${(molecule.properties.qed - originalProperties.qed).toFixed(2)}`}
                                      color={Math.abs(molecule.properties.qed - targetProperties.qed.target) < Math.abs(originalProperties.qed - targetProperties.qed.target) ? 'success' : 'error'}
                                      sx={{ ml: 1 }}
                                    />
                                  )}
                                </Typography>
                              </Grid>
                              <Grid item xs={6}>
                                <Typography variant="body2" color="text.secondary">
                                  Synthetic Accessibility:
                                </Typography>
                                <Typography variant="body2">
                                  {molecule.properties.synthetic_accessibility}
                                  {originalProperties && (
                                    <Chip 
                                      size="small" 
                                      label={`${(molecule.properties.synthetic_accessibility - originalProperties.synthetic_accessibility).toFixed(2)}`}
                                      color={Math.abs(molecule.properties.synthetic_accessibility - targetProperties.synthetic_accessibility.target) < Math.abs(originalProperties.synthetic_accessibility - targetProperties.synthetic_accessibility.target) ? 'success' : 'error'}
                                      sx={{ ml: 1 }}
                                    />
                                  )}
                                </Typography>
                              </Grid>
                              <Grid item xs={6}>
                                <Typography variant="body2" color="text.secondary">
                                  Solubility:
                                </Typography>
                                <Typography variant="body2">
                                  {molecule.properties.solubility}
                                  {originalProperties && (
                                    <Chip 
                                      size="small" 
                                      label={molecule.properties.solubility === targetProperties.solubility.target ? "✓" : "✗"}
                                      color={molecule.properties.solubility === targetProperties.solubility.target ? 'success' : 'error'}
                                      sx={{ ml: 1 }}
                                    />
                                  )}
                                </Typography>
                              </Grid>
                            </Grid>
                            
                            <Box sx={{ mt: 2 }}>
                              <Chip 
                                label={`Toxicity: ${molecule.properties.toxicity_risk}`} 
                                size="small" 
                                color={molecule.properties.toxicity_risk === 'Low' ? 'success' : molecule.properties.toxicity_risk === 'Moderate' ? 'warning' : 'error'}
                                sx={{ mr: 1, mb: 1 }} 
                              />
                              <Chip 
                                label={`Bioavailability: ${molecule.properties.bioavailability}`} 
                                size="small" 
                                sx={{ mr: 1, mb: 1 }} 
                              />
                            </Box>
                          </CardContent>
                        </Grid>
                      </Grid>
                    </Card>
                  )
                ))}
              </>
            )}
          </Paper>
        </Grid>
      </Grid>
    </Container>
  );
}

export default OptimizePage;