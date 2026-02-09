import { initializeApp } from "firebase/app";
import { getAnalytics } from "firebase/analytics";
import { getAuth, GoogleAuthProvider } from "firebase/auth";
import { getFirestore } from "firebase/firestore";
import { getStorage } from "firebase/storage";

const firebaseConfig = {
  apiKey: "AIzaSyAxYxzXTgrHWDXLWve6ldkNb4b6bKZ8eec",
  authDomain: "phoenix-4d0e9.firebaseapp.com",
  projectId: "phoenix-4d0e9",
  storageBucket: "phoenix-4d0e9.firebasestorage.app",
  messagingSenderId: "780546102728",
  appId: "1:780546102728:web:27359fd199ffcac325dcfd",
  measurementId: "G-V1E3QYG7MG"
};

// Initialize Firebase
const app = initializeApp(firebaseConfig);
const analytics = getAnalytics(app);
const auth = getAuth(app);
const db = getFirestore(app);
const storage = getStorage(app);
const googleProvider = new GoogleAuthProvider();

export { app, analytics, auth, db, storage, googleProvider };
