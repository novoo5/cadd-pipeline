import type { Metadata } from "next";
import { Inter } from "next/font/google";
import "./globals.css";

const inter = Inter({ subsets: ["latin"] });

export const metadata: Metadata = {
  title: "CADD Pipeline-Novoo",
  description: "Compound screening pipeline: analogue generation, ADMET filtering, molecular docking, and retrosynthesis.",
};

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en" className="dark">
      <body className={`${inter.className} bg-[#0a0a0a] text-gray-200 antialiased min-h-screen flex flex-col`}>
        {/* Minimal Nav */}
        <nav className="border-b border-gray-800/50 bg-[#0a0a0a]/80 backdrop-blur-sm sticky top-0 z-50">
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <div className="flex items-center justify-between h-14">
              <div className="flex items-center gap-3">
                <div className="w-7 h-7 rounded-md bg-emerald-900/40 border border-emerald-800 flex items-center justify-center">
                  <span className="text-emerald-500 font-bold text-xs">C</span>
                </div>
                <span className="font-medium text-gray-200 tracking-tight text-sm">
                  CADD Pipeline
                </span>
              </div>
            </div>
          </div>
        </nav>

        {/* Content */}
        <main className="flex-1 max-w-7xl w-full mx-auto px-4 sm:px-6 lg:px-8 py-10">
          {children}
        </main>
      </body>
    </html>
  );
}
