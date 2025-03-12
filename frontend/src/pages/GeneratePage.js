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
} from '@mui/material';
import {
  Science as ScienceIcon,
  Info as InfoIcon,
  Download as DownloadIcon,
  Bookmark as BookmarkIcon,
} from '@mui/icons-material';
import drugDiscoveryApi from '../services/api';

function GeneratePage() {
  // State for form inputs
  const [targetProperties, setTargetProperties] = useState({
    molecular_weight: { min: 200, max: 500 },
    logP: { min: 0, max: 5 },
    num_h_donors: { max: 5 },
    num_h_acceptors: { max: 10 },
  });
  
  const [constraints, setConstraints] = useState({
    include_substructure: '',
    exclude_substructure: '',
    must_have_elements: [],
  });
  
  const [numMolecules, setNumMolecules] = useState(5);
  
  // State for API response
  const [generatedMolecules, setGeneratedMolecules] = useState([]);
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
      const response = await drugDiscoveryApi.generateMolecules(
        targetProperties,
        constraints,
        numMolecules
      );
      
      setGeneratedMolecules(response.molecules);
      setSuccess(true);
      window.scrollTo({ top: document.getElementById('results').offsetTop, behavior: 'smooth' });
    } catch (error) {
      setError('Error generating molecules. Please try again.');
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
  
  // Handle constraint changes
  const handleConstraintChange = (constraint, value) => {
    setConstraints((prev) => ({
      ...prev,
      [constraint]: value,
    }));
  };
  
  return (
    <Container maxWidth="lg">
      <Typography variant="h4" component="h1" gutterBottom sx={{ mb: 3 }}>
        <ScienceIcon sx={{ mr: 1, verticalAlign: 'middle' }} />
        Generate Novel Drug Candidates
      </Typography>
      
      <Typography variant="body1" paragraph>
        Specify target properties and constraints to generate novel drug candidates using our AI-powered platform.
      </Typography>
      
      <Grid container spacing={4}>
        {/* Form Section */}
        <Grid item xs={12} md={5}>
          <Paper elevation={3} sx={{ p: 3, borderRadius: 2 }}>
            <form onSubmit={handleSubmit}>
              <Typography variant="h6" gutterBottom>
                Target Properties
              </Typography>
              
              <Box sx={{ mb: 3 }}>
                <Typography gutterBottom>
                  Molecular Weight (Da)
                </Typography>
                <Grid container spacing={2}>
                  <Grid item xs={6}>
                    <TextField
                      label="Min"
                      type="number"
                      value={targetProperties.molecular_weight.min}
                      onChange={(e) => handlePropertyChange('molecular_weight', 'min', Number(e.target.value))}
                      fullWidth
                      size="small"
                    />
                  </Grid>
                  <Grid item xs={6}>
                    <TextField
                      label="Max"
                      type="number"
                      value={targetProperties.molecular_weight.max}
                      onChange={(e) => handlePropertyChange('molecular_weight', 'max', Number(e.target.value))}
                      fullWidth
                      size="small"
                    />
                  </Grid>
                </Grid>
              </Box>
              
              <Box sx={{ mb: 3 }}>
                <Typography gutterBottom>
                  LogP (Lipophilicity)
                </Typography>
                <Grid container spacing={2}>
                  <Grid item xs={6}>
                    <TextField
                      label="Min"
                      type="number"
                      value={targetProperties.logP.min}
                      onChange={(e) => handlePropertyChange('logP', 'min', Number(e.target.value))}
                      fullWidth
                      size="small"
                    />
                  </Grid>
                  <Grid item xs={6}>
                    <TextField
                      label="Max"
                      type="number"
                      value={targetProperties.logP.max}
                      onChange={(e) => handlePropertyChange('logP', 'max', Number(e.target.value))}
                      fullWidth
                      size="small"
                    />
                  </Grid>
                </Grid>
              </Box>
              
              <Box sx={{ mb: 3 }}>
                <Typography gutterBottom>
                  Max Hydrogen Bond Donors
                </Typography>
                <Slider
                  value={targetProperties.num_h_donors.max}
                  onChange={(e, value) => handlePropertyChange('num_h_donors', 'max', value)}
                  valueLabelDisplay="auto"
                  step={1}
                  marks
                  min={0}
                  max={10}
                />
              </Box>
              
              <Box sx={{ mb: 3 }}>
                <Typography gutterBottom>
                  Max Hydrogen Bond Acceptors
                </Typography>
                <Slider
                  value={targetProperties.num_h_acceptors.max}
                  onChange={(e, value) => handlePropertyChange('num_h_acceptors', 'max', value)}
                  valueLabelDisplay="auto"
                  step={1}
                  marks
                  min={0}
                  max={20}
                />
              </Box>
              
              <Divider sx={{ my: 3 }} />
              
              <Typography variant="h6" gutterBottom>
                Constraints (Optional)
              </Typography>
              
              <Box sx={{ mb: 3 }}>
                <TextField
                  label="Include Substructure (SMILES)"
                  value={constraints.include_substructure}
                  onChange={(e) => handleConstraintChange('include_substructure', e.target.value)}
                  fullWidth
                  size="small"
                  placeholder="e.g., c1ccccc1"
                  helperText="SMILES notation of substructure that must be present"
                />
              </Box>
              
              <Box sx={{ mb: 3 }}>
                <TextField
                  label="Exclude Substructure (SMILES)"
                  value={constraints.exclude_substructure}
                  onChange={(e) => handleConstraintChange('exclude_substructure', e.target.value)}
                  fullWidth
                  size="small"
                  placeholder="e.g., C(=O)Cl"
                  helperText="SMILES notation of substructure that must not be present"
                />
              </Box>
              
              <Divider sx={{ my: 3 }} />
              
              <Box sx={{ mb: 3 }}>
                <Typography gutterBottom>
                  Number of Molecules to Generate
                </Typography>
                <Slider
                  value={numMolecules}
                  onChange={(e, value) => setNumMolecules(value)}
                  valueLabelDisplay="auto"
                  step={1}
                  marks
                  min={1}
                  max={10}
                />
              </Box>
              
              <Button
                type="submit"
                variant="contained"
                color="primary"
                fullWidth
                size="large"
                disabled={loading}
                startIcon={loading ? <CircularProgress size={24} /> : <ScienceIcon />}
              >
                {loading ? 'Generating...' : 'Generate Molecules'}
              </Button>
              
              {error && (
                <Alert severity="error" sx={{ mt: 2 }}>
                  {error}
                </Alert>
              )}
              
              {success && (
                <Alert severity="success" sx={{ mt: 2 }}>
                  Successfully generated {generatedMolecules.length} molecules!
                </Alert>
              )}
            </form>
          </Paper>
        </Grid>
        
        {/* Results Section */}
        <Grid item xs={12} md={7} id="results">
          <Paper elevation={3} sx={{ p: 3, borderRadius: 2 }}>
            <Typography variant="h6" gutterBottom>
              Generated Molecules
            </Typography>
            
            {generatedMolecules.length === 0 ? (
              <Box sx={{ textAlign: 'center', py: 4 }}>
                <ScienceIcon sx={{ fontSize: 60, color: 'text.secondary', opacity: 0.3 }} />
                <Typography variant="body1" color="text.secondary" sx={{ mt: 2 }}>
                  No molecules generated yet. Adjust the parameters and click "Generate Molecules" to start.
                </Typography>
              </Box>
            ) : (
              <Grid container spacing={3}>
                {generatedMolecules.map((molecule, index) => (
                  <Grid item xs={12} key={index}>
                    <Card sx={{ display: 'flex', flexDirection: { xs: 'column', sm: 'row' } }}>
                      <CardMedia
                        component="img"
                        sx={{ width: { xs: '100%', sm: 200 }, height: { xs: 200, sm: 'auto' } }}
                        image={molecule.image_url || 'https://via.placeholder.com/200?text=Molecule'}
                        alt={`Molecule ${index + 1}`}
                      />
                      <CardContent sx={{ flex: '1 0 auto' }}>
                        <Typography variant="h6" gutterBottom>
                          Molecule {index + 1}
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
                          {molecule.smiles}
                        </Typography>
                        
                        <Grid container spacing={1}>
                          <Grid item xs={6} sm={4}>
                            <Typography variant="body2" color="text.secondary">
                              Molecular Weight:
                            </Typography>
                            <Typography variant="body2">
                              {molecule.properties.molecular_weight} Da
                            </Typography>
                          </Grid>
                          <Grid item xs={6} sm={4}>
                            <Typography variant="body2" color="text.secondary">
                              LogP:
                            </Typography>
                            <Typography variant="body2">
                              {molecule.properties.logP}
                            </Typography>
                          </Grid>
                          <Grid item xs={6} sm={4}>
                            <Typography variant="body2" color="text.secondary">
                              H-Donors:
                            </Typography>
                            <Typography variant="body2">
                              {molecule.properties.num_h_donors}
                            </Typography>
                          </Grid>
                          <Grid item xs={6} sm={4}>
                            <Typography variant="body2" color="text.secondary">
                              H-Acceptors:
                            </Typography>
                            <Typography variant="body2">
                              {molecule.properties.num_h_acceptors}
                            </Typography>
                          </Grid>
                          <Grid item xs={6} sm={4}>
                            <Typography variant="body2" color="text.secondary">
                              QED:
                            </Typography>
                            <Typography variant="body2">
                              {molecule.properties.qed}
                            </Typography>
                          </Grid>
                          <Grid item xs={6} sm={4}>
                            <Typography variant="body2" color="text.secondary">
                              Druglikeness:
                            </Typography>
                            <Typography variant="body2">
                              {molecule.properties.druglikeness}
                            </Typography>
                          </Grid>
                        </Grid>
                        
                        <Box sx={{ mt: 2 }}>
                          <Chip 
                            label={`Solubility: ${molecule.properties.solubility}`} 
                            size="small" 
                            sx={{ mr: 1, mb: 1 }} 
                          />
                          <Chip 
                            label={`Bioavailability: ${molecule.properties.bioavailability}`} 
                            size="small" 
                            sx={{ mr: 1, mb: 1 }} 
                          />
                          <Chip 
                            label={`Toxicity: ${molecule.properties.toxicity_risk}`} 
                            size="small" 
                            color={molecule.properties.toxicity_risk === 'Low' ? 'success' : molecule.properties.toxicity_risk === 'Moderate' ? 'warning' : 'error'}
                            sx={{ mr: 1, mb: 1 }} 
                          />
                        </Box>
                      </CardContent>
                    </Card>
                  </Grid>
                ))}
              </Grid>
            )}
          </Paper>
        </Grid>
      </Grid>
    </Container>
  );
}

export default GeneratePage; 