import React from 'react';
import { Link as RouterLink } from 'react-router-dom';
import {
  Box,
  Button,
  Card,
  CardContent,
  CardMedia,
  Container,
  Grid,
  Typography,
  Paper,
} from '@mui/material';
import {
  Science as ScienceIcon,
  Tune as TuneIcon,
  Search as SearchIcon,
  Speed as SpeedIcon,
  Biotech as BiotechIcon,
  BarChart as BarChartIcon,
} from '@mui/icons-material';

const features = [
  {
    title: 'Generate Novel Molecules',
    description: 'Create new drug candidates based on desired properties using generative AI.',
    icon: <ScienceIcon fontSize="large" color="primary" />,
    path: '/generate',
  },
  {
    title: 'Optimize Existing Molecules',
    description: 'Improve existing drug candidates to enhance specific properties.',
    icon: <TuneIcon fontSize="large" color="primary" />,
    path: '/optimize',
  },
  {
    title: 'Search Similar Compounds',
    description: 'Find molecules similar to a query structure in our database.',
    icon: <SearchIcon fontSize="large" color="primary" />,
    path: '/search',
  },
];

const benefits = [
  {
    title: 'Accelerate Discovery',
    description: 'Reduce time to identify promising drug candidates from years to days.',
    icon: <SpeedIcon fontSize="large" color="secondary" />,
  },
  {
    title: 'Explore Chemical Space',
    description: 'Access a vast chemical space beyond traditional methods.',
    icon: <BiotechIcon fontSize="large" color="secondary" />,
  },
  {
    title: 'Data-Driven Decisions',
    description: 'Make informed decisions based on predicted properties and AI insights.',
    icon: <BarChartIcon fontSize="large" color="secondary" />,
  },
];

function HomePage() {
  return (
    <Container maxWidth="lg">
      {/* Hero Section */}
      <Paper
        sx={{
          position: 'relative',
          backgroundColor: 'grey.800',
          color: '#fff',
          mb: 4,
          backgroundSize: 'cover',
          backgroundRepeat: 'no-repeat',
          backgroundPosition: 'center',
          backgroundImage: 'url(https://source.unsplash.com/random?science,molecule)',
          borderRadius: 2,
          overflow: 'hidden',
        }}
      >
        <Box
          sx={{
            position: 'absolute',
            top: 0,
            bottom: 0,
            right: 0,
            left: 0,
            backgroundColor: 'rgba(0,0,0,.6)',
          }}
        />
        <Grid container>
          <Grid item md={6}>
            <Box
              sx={{
                position: 'relative',
                p: { xs: 3, md: 6 },
                pr: { md: 0 },
              }}
            >
              <Typography component="h1" variant="h3" color="inherit" gutterBottom>
                Accelerating Drug Discovery with AI
              </Typography>
              <Typography variant="h5" color="inherit" paragraph>
                Leverage the power of generative AI to design, optimize, and discover novel drug candidates faster than ever before.
              </Typography>
              <Button
                variant="contained"
                size="large"
                component={RouterLink}
                to="/generate"
                sx={{ mt: 2 }}
              >
                Get Started
              </Button>
            </Box>
          </Grid>
        </Grid>
      </Paper>

      {/* Features Section */}
      <Typography variant="h4" component="h2" gutterBottom align="center" sx={{ mt: 8, mb: 4 }}>
        Key Features
      </Typography>
      <Grid container spacing={4}>
        {features.map((feature, index) => (
          <Grid item key={index} xs={12} sm={6} md={4}>
            <Card
              sx={{
                height: '100%',
                display: 'flex',
                flexDirection: 'column',
                transition: 'transform 0.3s ease-in-out, box-shadow 0.3s ease-in-out',
                '&:hover': {
                  transform: 'translateY(-8px)',
                  boxShadow: '0 12px 20px rgba(0, 0, 0, 0.2)',
                },
              }}
            >
              <CardContent sx={{ flexGrow: 1, textAlign: 'center', pt: 4 }}>
                <Box sx={{ mb: 2 }}>{feature.icon}</Box>
                <Typography gutterBottom variant="h5" component="h3">
                  {feature.title}
                </Typography>
                <Typography variant="body1" color="text.secondary" sx={{ mb: 3 }}>
                  {feature.description}
                </Typography>
                <Button
                  variant="outlined"
                  color="primary"
                  component={RouterLink}
                  to={feature.path}
                >
                  Learn More
                </Button>
              </CardContent>
            </Card>
          </Grid>
        ))}
      </Grid>

      {/* Benefits Section */}
      <Box sx={{ bgcolor: 'background.paper', py: 8, mt: 8, borderRadius: 2 }}>
        <Container maxWidth="lg">
          <Typography variant="h4" component="h2" gutterBottom align="center" sx={{ mb: 6 }}>
            Benefits
          </Typography>
          <Grid container spacing={4}>
            {benefits.map((benefit, index) => (
              <Grid item key={index} xs={12} md={4}>
                <Box sx={{ textAlign: 'center' }}>
                  <Box sx={{ mb: 2 }}>{benefit.icon}</Box>
                  <Typography variant="h5" component="h3" gutterBottom>
                    {benefit.title}
                  </Typography>
                  <Typography variant="body1" color="text.secondary">
                    {benefit.description}
                  </Typography>
                </Box>
              </Grid>
            ))}
          </Grid>
        </Container>
      </Box>

      {/* Call to Action */}
      <Box
        sx={{
          bgcolor: 'primary.main',
          color: 'white',
          py: 6,
          mt: 8,
          borderRadius: 2,
          textAlign: 'center',
        }}
      >
        <Typography variant="h4" component="h2" gutterBottom>
          Ready to Revolutionize Your Drug Discovery Process?
        </Typography>
        <Typography variant="body1" paragraph sx={{ maxWidth: 600, mx: 'auto', mb: 4 }}>
          Start using our AI-powered platform today to accelerate your research and discover promising drug candidates.
        </Typography>
        <Button
          variant="contained"
          color="secondary"
          size="large"
          component={RouterLink}
          to="/generate"
        >
          Get Started Now
        </Button>
      </Box>
    </Container>
  );
}

export default HomePage; 