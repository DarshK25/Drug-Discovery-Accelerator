import React from 'react';
import {
  Box,
  Container,
  Typography,
  Paper,
  Grid,
  Card,
  CardContent,
  CardMedia,
  Divider,
  List,
  ListItem,
  ListItemIcon,
  ListItemText,
} from '@mui/material';
import {
  Science as ScienceIcon,
  Speed as SpeedIcon,
  Biotech as BiotechIcon,
  Psychology as PsychologyIcon,
  DataObject as DataObjectIcon,
  Lightbulb as LightbulbIcon,
} from '@mui/icons-material';

function AboutPage() {
  return (
    <Container maxWidth="lg">
      <Typography variant="h4" component="h1" gutterBottom sx={{ mb: 3 }}>
        About Drug Discovery Assistant
      </Typography>
      
      <Paper elevation={3} sx={{ p: 4, mb: 4, borderRadius: 2 }}>
        <Typography variant="h5" gutterBottom>
          Our Mission
        </Typography>
        <Typography variant="body1" paragraph>
          The Drug Discovery Assistant aims to revolutionize the drug discovery process by leveraging the power of generative AI. 
          Our platform enables researchers to efficiently identify, design, and optimize potential drug candidates, 
          reducing the time and cost associated with traditional drug discovery methods.
        </Typography>
        <Typography variant="body1" paragraph>
          Traditional drug discovery is time-consuming, expensive, and often inefficient, with a high rate of failure in clinical trials. 
          By combining advanced AI techniques with cheminformatics, we provide tools that can accelerate the early stages of drug discovery, 
          allowing researchers to focus on the most promising candidates.
        </Typography>
      </Paper>
      
      <Grid container spacing={4} sx={{ mb: 4 }}>
        <Grid item xs={12} md={6}>
          <Card sx={{ height: '100%' }}>
            <CardMedia
              component="img"
              height="240"
              image="https://source.unsplash.com/random?science,laboratory"
              alt="Laboratory"
            />
            <CardContent>
              <Typography variant="h5" gutterBottom>
                The Challenge
              </Typography>
              <Typography variant="body1" color="text.secondary">
                Drug discovery traditionally takes 10-15 years and costs billions of dollars, with a high failure rate. 
                The vast chemical space makes it challenging to identify promising candidates efficiently. 
                Our platform addresses these challenges by using AI to explore chemical space more effectively and predict properties accurately.
              </Typography>
            </CardContent>
          </Card>
        </Grid>
        <Grid item xs={12} md={6}>
          <Card sx={{ height: '100%' }}>
            <CardMedia
              component="img"
              height="240"
              image="https://source.unsplash.com/random?technology,ai"
              alt="AI Technology"
            />
            <CardContent>
              <Typography variant="h5" gutterBottom>
                Our Solution
              </Typography>
              <Typography variant="body1" color="text.secondary">
                The Drug Discovery Assistant combines generative AI models, molecular property prediction, and similarity search 
                to provide a comprehensive platform for drug discovery. Our tools enable researchers to generate novel compounds, 
                optimize existing ones, and identify promising candidates with desired properties.
              </Typography>
            </CardContent>
          </Card>
        </Grid>
      </Grid>
      
      <Typography variant="h5" gutterBottom sx={{ mb: 3 }}>
        Key Technologies
      </Typography>
      
      <Grid container spacing={3} sx={{ mb: 4 }}>
        <Grid item xs={12} sm={6} md={4}>
          <Card sx={{ height: '100%' }}>
            <CardContent>
              <Box sx={{ display: 'flex', alignItems: 'center', mb: 2 }}>
                <ScienceIcon color="primary" sx={{ fontSize: 40, mr: 2 }} />
                <Typography variant="h6">
                  RDKit
                </Typography>
              </Box>
              <Typography variant="body2" color="text.secondary">
                An open-source cheminformatics toolkit that provides comprehensive functionality for working with chemical structures, 
                including molecule manipulation, property calculation, and visualization.
              </Typography>
            </CardContent>
          </Card>
        </Grid>
        <Grid item xs={12} sm={6} md={4}>
          <Card sx={{ height: '100%' }}>
            <CardContent>
              <Box sx={{ display: 'flex', alignItems: 'center', mb: 2 }}>
                <PsychologyIcon color="primary" sx={{ fontSize: 40, mr: 2 }} />
                <Typography variant="h6">
                  PyTorch
                </Typography>
              </Box>
              <Typography variant="body2" color="text.secondary">
                A deep learning framework that powers our generative models, enabling the creation of novel molecular structures 
                and the prediction of molecular properties with high accuracy.
              </Typography>
            </CardContent>
          </Card>
        </Grid>
        <Grid item xs={12} sm={6} md={4}>
          <Card sx={{ height: '100%' }}>
            <CardContent>
              <Box sx={{ display: 'flex', alignItems: 'center', mb: 2 }}>
                <DataObjectIcon color="primary" sx={{ fontSize: 40, mr: 2 }} />
                <Typography variant="h6">
                  Transformers
                </Typography>
              </Box>
              <Typography variant="body2" color="text.secondary">
                State-of-the-art natural language processing models adapted for molecular generation, 
                allowing our platform to learn from vast datasets of chemical structures and generate novel compounds.
              </Typography>
            </CardContent>
          </Card>
        </Grid>
      </Grid>
      
      <Divider sx={{ my: 4 }} />
      
      <Typography variant="h5" gutterBottom sx={{ mb: 3 }}>
        Benefits for Drug Discovery
      </Typography>
      
      <List>
        <ListItem>
          <ListItemIcon>
            <SpeedIcon color="secondary" />
          </ListItemIcon>
          <ListItemText 
            primary="Accelerated Discovery" 
            secondary="Reduce the time to identify promising drug candidates from years to months or even weeks." 
          />
        </ListItem>
        <ListItem>
          <ListItemIcon>
            <BiotechIcon color="secondary" />
          </ListItemIcon>
          <ListItemText 
            primary="Expanded Chemical Space" 
            secondary="Explore a vast chemical space beyond what traditional methods can access, increasing the chances of finding novel compounds." 
          />
        </ListItem>
        <ListItem>
          <ListItemIcon>
            <LightbulbIcon color="secondary" />
          </ListItemIcon>
          <ListItemText 
            primary="Intelligent Design" 
            secondary="Use AI-driven insights to design molecules with specific properties, optimizing for efficacy, safety, and synthesizability." 
          />
        </ListItem>
      </List>
      
      <Box sx={{ bgcolor: 'background.paper', p: 4, mt: 4, borderRadius: 2, textAlign: 'center' }}>
        <Typography variant="h5" gutterBottom>
          Join Us in Revolutionizing Drug Discovery
        </Typography>
        <Typography variant="body1">
          We're committed to advancing the field of drug discovery through innovative AI technologies. 
          Our platform is continuously evolving, with new features and improvements being added regularly.
        </Typography>
      </Box>
    </Container>
  );
}

export default AboutPage;