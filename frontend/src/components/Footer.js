import React from 'react';
import { Box, Container, Typography, Link, Grid } from '@mui/material';
import { GitHub as GitHubIcon, LinkedIn as LinkedInIcon, Twitter as TwitterIcon } from '@mui/icons-material';

function Footer() {
  return (
    <Box
      component="footer"
      sx={{
        py: 3,
        px: 2,
        mt: 'auto',
        backgroundColor: (theme) => theme.palette.grey[900],
        color: 'white',
      }}
    >
      <Container maxWidth="lg">
        <Grid container spacing={3}>
          <Grid item xs={12} sm={4}>
            <Typography variant="h6" gutterBottom>
              Drug Discovery Assistant
            </Typography>
            <Typography variant="body2" color="text.secondary" sx={{ color: 'rgba(255, 255, 255, 0.7)' }}>
              Accelerating drug discovery with generative AI
            </Typography>
          </Grid>
          <Grid item xs={12} sm={4}>
            <Typography variant="h6" gutterBottom>
              Resources
            </Typography>
            <Typography variant="body2">
              <Link href="https://www.rdkit.org/" target="_blank" rel="noopener" sx={{ color: 'rgba(255, 255, 255, 0.7)' }}>
                RDKit
              </Link>
            </Typography>
            <Typography variant="body2">
              <Link href="https://pytorch.org/" target="_blank" rel="noopener" sx={{ color: 'rgba(255, 255, 255, 0.7)' }}>
                PyTorch
              </Link>
            </Typography>
            <Typography variant="body2">
              <Link href="https://huggingface.co/" target="_blank" rel="noopener" sx={{ color: 'rgba(255, 255, 255, 0.7)' }}>
                Hugging Face
              </Link>
            </Typography>
          </Grid>
          <Grid item xs={12} sm={4}>
            <Typography variant="h6" gutterBottom>
              Connect
            </Typography>
            <Box sx={{ display: 'flex', gap: 2 }}>
              <Link href="https://github.com/" target="_blank" rel="noopener" sx={{ color: 'rgba(255, 255, 255, 0.7)' }}>
                <GitHubIcon />
              </Link>
              <Link href="https://linkedin.com/" target="_blank" rel="noopener" sx={{ color: 'rgba(255, 255, 255, 0.7)' }}>
                <LinkedInIcon />
              </Link>
              <Link href="https://twitter.com/" target="_blank" rel="noopener" sx={{ color: 'rgba(255, 255, 255, 0.7)' }}>
                <TwitterIcon />
              </Link>
            </Box>
          </Grid>
        </Grid>
        <Box mt={3}>
          <Typography variant="body2" color="text.secondary" align="center" sx={{ color: 'rgba(255, 255, 255, 0.7)' }}>
            {'Copyright © '}
            <Link color="inherit" href="/">
              Drug Discovery Assistant
            </Link>{' '}
            {new Date().getFullYear()}
            {'.'}
          </Typography>
        </Box>
      </Container>
    </Box>
  );
}

export default Footer; 