import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  // Pass the IP address without the port
  allowedDevOrigins: [
    "192.168.1.4",
    "localhost",
    "127.0.0.1"
  ],

  typescript: {
    ignoreBuildErrors: true,
  }
};

export default nextConfig;
