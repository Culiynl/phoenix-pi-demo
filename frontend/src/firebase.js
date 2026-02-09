import { initializeApp } from "firebase/app";
import { getAnalytics } from "firebase/analytics";
import { getAuth, GoogleAuthProvider } from "firebase/auth";
import { getFirestore } from "firebase/firestore";
import { getStorage } from "firebase/storage";

const firebaseConfig = {
  apiKey: "AIzaSyBZywnI16E3EjBaTlSLrsi1abD5Rmxu_5o",
  authDomain: "pheonix-web-53172.firebaseapp.com",
  projectId: "pheonix-web-53172",
  storageBucket: "pheonix-web-53172.firebasestorage.app",
  messagingSenderId: "236113511528",
  appId: "1:236113511528:web:a98f2de9c9ca11871c67dd",
  measurementId: "G-9CJQ2R7JK2"
};

// Initialize Firebase
const app = initializeApp(firebaseConfig);
const analytics = getAnalytics(app);
const auth = getAuth(app);
const db = getFirestore(app);
const storage = getStorage(app);
const googleProvider = new GoogleAuthProvider();

export { app, analytics, auth, db, storage, googleProvider };
